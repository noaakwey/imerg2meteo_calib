"""
apply_qm_to_rasters.py
======================
Применяет QM-калибровку к архивам квартальных растров IMERG или ERA5-Land,
группирует кварталы по годам, сохраняет годовые откалиброванные GeoTIFF.

Поддерживаемые датасеты
-----------------------
  imerg
      Источник: NASA GPM IMERG V07 (GEE export)
      Файлы   : IMERG_V07_P30min_mmh_YYYY_QN_permanent.tif
      Шаг     : 30 мин (~4320 бэндов/квартал)
      Единицы : мм/ч (интенсивность за 30-мин слот)
      Бэнды   : 'P_YYYYMMDD_HHmm'

  era5land
      Источник: ECMWF ERA5-Land Hourly (GEE export)
      Файлы   : ERA5Land_P1h_mm_YYYY_QN.tif
      Шаг     : 1 ч (~2160 бэндов/квартал)
      Единицы : мм/ч (total_precipitation_hourly × 1000)
      Бэнды   : 'P_YYYYMMDD_HHmm'

Алгоритм
---------
1. Распаковка ZIP во временный каталог, группировка файлов по годам.
2. Объединение бэндов квартальных файлов в хронологически упорядоченный
   годовой ряд.
3. KNN-сопоставление пикселей растра с ближайшей метеостанцией (cKDTree).
4. Для каждого 3-часового окна:
      a. Пересчёт интенсивности → накопленные мм (интенсивность × шаг_ч)
      b. Суммирование по слотам → 3ч сумма осадков (мм)
      c. QM-коррекция (Full-Distribution, 1000 квантилей + Volume Scaling)
         отдельно для каждой зоны Вороного (ближайшей станции)
      d. Пропорциональная дезагрегация откалиброванной 3ч суммы обратно
         по исходным субвременны́м слотам
      e. Обратный пересчёт накопленных мм → интенсивность мм/ч
5. Запись годового GeoTIFF (float32, LZW, tiled 256×256).

Использование
-------------
  python apply_qm_to_rasters.py --dataset imerg
  python apply_qm_to_rasters.py --dataset era5land
  python apply_qm_to_rasters.py --dataset imerg --year 2010
  python apply_qm_to_rasters.py --dataset imerg \\
      --zip  "path/to/archive.zip" \\
      --out  "path/to/output/" \\
      --calib "path/to/calib_imerg/"
"""

import os
import re
import sys
import glob
import zipfile
import tempfile
import argparse
import shutil
from collections import defaultdict

# Принудительно UTF-8 для stdout/stderr (Windows CP1251 не содержит ряд символов)
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
if hasattr(sys.stderr, 'reconfigure'):
    sys.stderr.reconfigure(encoding='utf-8', errors='replace')

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from tqdm import tqdm

try:
    import rioxarray as rxr
    import rasterio
except ImportError:
    print("ERROR: Требуются пакеты rioxarray и rasterio.")
    print("  pip install rioxarray rasterio")
    sys.exit(1)

# ------------------------------------------------------------------
# Добавляем src в sys.path, чтобы импортировать qm_calibration
# ------------------------------------------------------------------
_SRC = os.path.dirname(os.path.abspath(__file__))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from qm_calibration import get_season, apply_qm


# ==================================================================
# Конфигурация наборов данных
# ==================================================================

DATASET_CONFIG = {
    'imerg': {
        # Регулярное выражение для извлечения года и квартала из имени файла
        'file_re': re.compile(
            r'IMERG_V07_P30min_mmh_(\d{4})_(Q\d)_permanent\.tif',
            re.IGNORECASE
        ),
        # Шаг в минутах
        'step_min': 30,
        # Слотов в одном 3-часовом окне QM
        'slots_per_3h': 6,
        # Каталог с файлами *_calib.csv (относительно BASE_DIR)
        'calib_subdir': os.path.join('output', 'calib_imerg'),
        # Подкаталог вывода (относительно BASE_DIR)
        'out_subdir': os.path.join('output', 'imerg_rfactor_calib_v5_year_anchor'),
        # Шаблон имени выходного файла: {year}
        'out_fname_tpl': 'IMERG_V07_P30min_mmh_{year}_calib_qm.tif',
        # Тренировочный период QM
        'train_start': '2001-01-01',
        'train_end':   '2015-12-31',
    },
    'era5land': {
        'file_re': re.compile(
            r'ERA5Land_P1h_mm_(\d{4})_(Q\d)\.tif',
            re.IGNORECASE
        ),
        'step_min': 60,
        'slots_per_3h': 3,
        'calib_subdir': os.path.join('output', 'calib_era5land'),
        'out_subdir': os.path.join('output', 'era5land_rfactor_calib'),
        'out_fname_tpl': 'ERA5Land_P1h_mm_{year}_calib_qm.tif',
        'train_start': '2001-01-01',
        'train_end':   '2015-12-31',
    },
}


# ==================================================================
# 1. Загрузка метаданных станций
# ==================================================================

def _expected_3h_steps(year: int) -> int:
    """Expected number of 3-hour synoptic rows for a full UTC year."""
    leap = (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0)
    return (366 if leap else 365) * 8


def load_station_metadata(meteo_dir: str) -> pd.DataFrame:
    """Reads WMO index and station coordinates from CSV headers (X=lat, Y=lon)."""
    csv_files = glob.glob(os.path.join(meteo_dir, '*.csv'))
    meta = []
    for f in csv_files:
        try:
            df = pd.read_csv(f, sep=';', encoding='cp866', nrows=1)
            meta.append({
                'wmo_index': int(df['Index'].iloc[0]),
                'name':      str(df['StationName'].iloc[0]),
                'lon':       float(df['Y'].iloc[0]),
                'lat':       float(df['X'].iloc[0]),
            })
        except Exception:
            pass
    if not meta:
        raise RuntimeError(f'Не найдено ни одной метеостанции в {meteo_dir}')
    return pd.DataFrame(meta)


# ==================================================================
# 2. Предварительный расчёт QM-моделей из calib CSV
# ==================================================================

def _fit_empirical_qm(p_sat: np.ndarray,
                      p_st: np.ndarray,
                      q_levels: np.ndarray,
                      min_samples: int = 30):
    """Returns (q_sat, q_st, slope, p_th) or None if data are insufficient."""
    if len(p_sat) < min_samples or len(p_st) < min_samples:
        return None

    q_sat = np.quantile(p_sat, q_levels)
    q_st = np.quantile(p_st, q_levels)

    if (q_sat[-1] - q_sat[-5]) > 0:
        slope = (q_st[-1] - q_st[-5]) / (q_sat[-1] - q_sat[-5])
    else:
        slope = 1.0

    return q_sat, q_st, slope, 0.0


def _calc_pbias(sim: np.ndarray, obs: np.ndarray) -> float:
    mask = np.isfinite(sim) & np.isfinite(obs)
    sim = sim[mask]
    obs = obs[mask]
    if len(obs) == 0:
        return np.nan
    denom = float(np.sum(obs))
    if denom == 0:
        return np.nan
    return 100.0 * float(np.sum(sim - obs)) / denom


def _calc_kge(sim: np.ndarray, obs: np.ndarray) -> float:
    mask = np.isfinite(sim) & np.isfinite(obs)
    sim = sim[mask]
    obs = obs[mask]
    if len(obs) < 2 or np.std(obs) == 0:
        return np.nan
    r = np.corrcoef(obs, sim)[0, 1]
    alpha = np.std(sim) / np.std(obs)
    beta = np.mean(sim) / np.mean(obs)
    return 1.0 - np.sqrt((r - 1) ** 2 + (alpha - 1) ** 2 + (beta - 1) ** 2)


def _score_alpha_candidate(df_tune: pd.DataFrame, candidate: np.ndarray) -> float:
    """Composite score: prefer daily/monthly dynamics with mild bias penalties."""
    tmp = df_tune[['datetime', 'P_station_mm']].copy()
    tmp['P_candidate_mm'] = candidate

    tmp['date'] = tmp['datetime'].dt.floor('D')
    daily = tmp.groupby('date', as_index=False)[['P_station_mm', 'P_candidate_mm']].sum()
    daily_kge = _calc_kge(daily['P_candidate_mm'].values, daily['P_station_mm'].values)
    daily_pbias = _calc_pbias(daily['P_candidate_mm'].values, daily['P_station_mm'].values)

    tmp['ym'] = tmp['datetime'].dt.to_period('M')
    monthly = tmp.groupby('ym', as_index=False)[['P_station_mm', 'P_candidate_mm']].sum()
    monthly_kge = _calc_kge(monthly['P_candidate_mm'].values, monthly['P_station_mm'].values)
    monthly_pbias = _calc_pbias(monthly['P_candidate_mm'].values, monthly['P_station_mm'].values)

    tmp['year'] = tmp['datetime'].dt.year
    annual = tmp.groupby('year', as_index=False)[['P_station_mm', 'P_candidate_mm']].sum()
    annual_pbias = _calc_pbias(annual['P_candidate_mm'].values, annual['P_station_mm'].values)

    dk = -2.0 if pd.isna(daily_kge) else daily_kge
    mk = -2.0 if pd.isna(monthly_kge) else monthly_kge
    dp = 200.0 if pd.isna(daily_pbias) else abs(daily_pbias)
    mp = 200.0 if pd.isna(monthly_pbias) else abs(monthly_pbias)
    ap = 200.0 if pd.isna(annual_pbias) else abs(annual_pbias)

    return (0.65 * mk + 0.45 * dk - 0.003 * mp - 0.002 * dp - 0.001 * ap)


def _fit_station_seasonal_models(train_df: pd.DataFrame,
                                 q_levels: np.ndarray,
                                 min_samples: int = 25) -> dict:
    """Fit seasonal QM models for one station: season -> (q_sat, q_st, slope, p_th)."""
    out = {}
    for season in ('DJF', 'MAM', 'JJA', 'SON'):
        sd = train_df[train_df['season'] == season]
        params = _fit_empirical_qm(
            sd['P_sat_mm'].values,
            sd['P_station_mm'].values,
            q_levels=q_levels,
            min_samples=min_samples,
        )
        if params is not None:
            out[season] = params
    return out


def _apply_station_seasonal_qm(df: pd.DataFrame, season_models: dict) -> np.ndarray:
    """Apply station seasonal models to a station dataframe and return corrected array."""
    corrected = df['P_sat_mm'].to_numpy(dtype=float, copy=True)
    for season, params in season_models.items():
        mask = (df['season'] == season).to_numpy()
        if not np.any(mask):
            continue
        q_sat, q_st, slope, p_th = params
        mapped = apply_qm(df.loc[mask, 'P_sat_mm'], q_sat, q_st, slope, p_th)
        corrected[mask] = np.maximum(mapped, 0.0)
    return corrected


def _select_station_blend_alpha(train_df: pd.DataFrame,
                                q_levels: np.ndarray,
                                default_alpha: float = 0.2) -> float:
    """
    Estimate station blend alpha for raster application:
      corrected = raw + alpha * (qm - raw)
    Uses blocked split inside train period.
    """
    if len(train_df) < 2500:
        return float(default_alpha)

    split_date = train_df['datetime'].quantile(0.75)
    subtrain = train_df[train_df['datetime'] <= split_date].copy()
    tune = train_df[train_df['datetime'] > split_date].copy()
    if len(subtrain) < 1500 or len(tune) < 500:
        return float(default_alpha)

    models_sub = _fit_station_seasonal_models(subtrain, q_levels=q_levels, min_samples=25)
    if len(models_sub) < 2:
        return float(default_alpha)

    qm_tune_full = _apply_station_seasonal_qm(train_df, models_sub)
    qm_tune = qm_tune_full[train_df['datetime'] > split_date]
    raw_tune = tune['P_sat_mm'].to_numpy(dtype=float)

    alpha_grid = [0.0, 0.1, 0.2, 0.3, 0.4]
    best_alpha = float(default_alpha)
    best_score = -1e18
    for alpha in alpha_grid:
        cand = raw_tune + alpha * (qm_tune - raw_tune)
        score = _score_alpha_candidate(tune, cand)
        if score > best_score:
            best_score = score
            best_alpha = float(alpha)

    return best_alpha


def precalculate_qm_models(calib_dir: str,
                            train_start: str,
                            train_end: str) -> dict:
    """
    Builds station models from *_calib.csv:
      1) seasonal 3h QM models (main correction),
      2) daily QM models (post-disaggregation volume control),
      3) annual ratio statistics (sanity guard).
    """
    seasonal_models = {s: {} for s in ('DJF', 'MAM', 'JJA', 'SON')}
    daily_models = {}
    annual_stats = {}
    annual_year_ratio = {}
    annual_transfer = {}
    blend_alpha = {}
    calib_files = glob.glob(os.path.join(calib_dir, '*_calib.csv'))

    if not calib_files:
        raise RuntimeError(f'No *_calib.csv files found in {calib_dir}')

    num_q = 1000
    q_levels = np.linspace(0.0, 1.0, num_q + 2)[1:-1]

    for f in tqdm(calib_files, desc='Training QM models'):
        try:
            df = pd.read_csv(f)
        except Exception as e:
            print(f'  Skip {f}: {e}')
            continue

        if 'wmo_index' not in df.columns:
            continue

        wmo = int(df['wmo_index'].iloc[0])
        df['datetime'] = pd.to_datetime(df['datetime'])
        if 'season' not in df.columns:
            df['season'] = df['datetime'].dt.month.apply(get_season)

        train_mask = (df['datetime'] >= train_start) & (df['datetime'] <= train_end)
        train_df = df[train_mask].copy()
        if train_df.empty:
            continue

        # 1) Seasonal 3h QM
        station_models = {}
        for season in seasonal_models:
            sd = train_df[train_df['season'] == season]
            params = _fit_empirical_qm(
                sd['P_sat_mm'].values,
                sd['P_station_mm'].values,
                q_levels=q_levels,
                min_samples=30,
            )
            if params is not None:
                seasonal_models[season][wmo] = params
                station_models[season] = params

        # Station-specific blend alpha (soft correction strength)
        if station_models:
            blend_alpha[wmo] = _select_station_blend_alpha(
                train_df,
                q_levels=q_levels,
                default_alpha=0.2,
            )

        # 2) Daily QM (sum over UTC day)
        daily_df = (
            train_df
            .assign(date=train_df['datetime'].dt.floor('D'))
            .groupby('date', as_index=False)[['P_sat_mm', 'P_station_mm']]
            .sum()
        )
        daily_params = _fit_empirical_qm(
            daily_df['P_sat_mm'].values,
            daily_df['P_station_mm'].values,
            q_levels=q_levels,
            min_samples=20,
        )
        if daily_params is not None:
            daily_models[wmo] = daily_params

        # 3) Annual ratio stats for final sanity guard
        annual_df = (
            train_df
            .assign(year=train_df['datetime'].dt.year)
            .groupby('year', as_index=False)[['P_sat_mm', 'P_station_mm']]
            .sum()
        )
        annual_df = annual_df[annual_df['P_sat_mm'] > 0]

        if len(annual_df) >= 3:
            ratio = (annual_df['P_station_mm'] / annual_df['P_sat_mm']).values
            ratio = np.clip(ratio, 0.2, 5.0)
            annual_stats[wmo] = {
                'ratio_p10': float(np.percentile(ratio, 10)),
                'ratio_p50': float(np.percentile(ratio, 50)),
                'ratio_p90': float(np.percentile(ratio, 90)),
                'station_p90': float(np.percentile(annual_df['P_station_mm'].values, 90)),
            }

        # 4) Year-specific annual anchor ratios (full years only, all available data)
        annual_full = (
            df
            .assign(year=df['datetime'].dt.year, n_rows=1)
            .groupby('year', as_index=False)[['P_sat_mm', 'P_station_mm', 'n_rows']]
            .sum()
        )
        ratio_by_year = {}
        for _, yr_row in annual_full.iterrows():
            y = int(yr_row['year'])
            n_rows = int(yr_row['n_rows'])
            expected = _expected_3h_steps(y)
            if expected <= 0 or n_rows < int(0.99 * expected):
                continue

            sat_sum = float(yr_row['P_sat_mm'])
            st_sum = float(yr_row['P_station_mm'])
            if sat_sum <= 0 or st_sum < 0:
                continue

            ratio_by_year[y] = float(np.clip(st_sum / max(sat_sum, 1e-6), 0.2, 5.0))

        if ratio_by_year:
            annual_year_ratio[wmo] = ratio_by_year

        # 5) Annual raw->station transfer model (targets interannual extremes)
        transfer_model = _fit_annual_transfer_model(annual_df)
        if transfer_model is not None:
            annual_transfer[wmo] = transfer_model

    total_seasonal = sum(len(v) for v in seasonal_models.values())
    print(f'  Done. Seasonal pairs (season x station): {total_seasonal}')
    print(f'  Daily models: {len(daily_models)} stations')
    print(f'  Annual guard stats: {len(annual_stats)} stations')
    print(f'  Year-anchor annual ratios: {len(annual_year_ratio)} stations')
    print(f'  Annual transfer models: {len(annual_transfer)} stations')
    if annual_year_ratio:
        n_years = np.array([len(v) for v in annual_year_ratio.values()], dtype=np.float64)
        print(f'  Year-anchor coverage: median={np.median(n_years):.1f} years/station, '
              f'min={np.min(n_years):.0f}, max={np.max(n_years):.0f}')
    if blend_alpha:
        alphas = np.array(list(blend_alpha.values()), dtype=float)
        print(f'  Blend alpha stats: median={np.median(alphas):.2f}, mean={np.mean(alphas):.2f}, '
              f'min={np.min(alphas):.2f}, max={np.max(alphas):.2f}')

    return {
        'seasonal': seasonal_models,
        'daily': daily_models,
        'annual': annual_stats,
        'annual_year_ratio': annual_year_ratio,
        'annual_transfer': annual_transfer,
        'blend_alpha': blend_alpha,
    }


# ==================================================================
# 3. Разбор имён бэндов → datetime
# ==================================================================

_BAND_PAT = re.compile(r'P_(\d{8})_(\d{4})')

def parse_band_datetimes(long_name_tuple) -> list:
    """
    'N_P_YYYYMMDD_HHmm' → list[pd.Timestamp | None].
    """
    dts = []
    for name in long_name_tuple:
        m = _BAND_PAT.search(str(name))
        if m:
            dts.append(pd.to_datetime(m.group(1) + m.group(2), format='%Y%m%d%H%M'))
        else:
            dts.append(None)
    return dts


# ==================================================================
# 4. QM-коррекция одного 3ч окна (2D-массивы)
# ==================================================================

def apply_qm_to_3h_window(raw_3h_mm: np.ndarray,
                           station_masks: dict,
                           season: str,
                           seasonal_models: dict,
                           blend_alpha: dict | None = None) -> np.ndarray:
    """
    Applies seasonal 3h QM correction to a 2D field of 3-hour accumulations (mm).
    Pixels without a station model remain unchanged.
    """
    calib = np.copy(raw_3h_mm).astype(np.float64)
    season_models = seasonal_models.get(season, {})

    for wmo, mask in station_masks.items():
        vals = raw_3h_mm[mask]
        if len(vals) == 0:
            continue

        params = season_models.get(wmo)
        if params is None:
            continue

        q_sat, q_st, slope, p_th = params
        corrected = np.maximum(apply_qm(vals, q_sat, q_st, slope, p_th), 0.0)

        # Soft station-specific blending keeps results close to synoptic signal.
        alpha = 0.2 if blend_alpha is None else float(blend_alpha.get(int(wmo), 0.2))
        alpha = float(np.clip(alpha, 0.0, 0.6))
        blended = vals + alpha * (corrected - vals)
        calib[mask] = np.maximum(blended, 0.0)

    return calib


def _build_day_slices(dts_arr: list) -> list:
    """Builds contiguous [i0, i1) slices for each UTC day in sorted dts_arr."""
    if not dts_arr:
        return []

    day_slices = []
    i0 = 0
    while i0 < len(dts_arr):
        d0 = dts_arr[i0].date()
        i1 = i0 + 1
        while i1 < len(dts_arr) and dts_arr[i1].date() == d0:
            i1 += 1
        day_slices.append((i0, i1))
        i0 = i1
    return day_slices


def apply_daily_qm_constraint(calib_mm: np.ndarray,
                              raw_mm: np.ndarray,
                              dts_arr: list,
                              station_masks: dict,
                              daily_models: dict) -> np.ndarray:
    """
    Daily post-correction: maps raw daily sum -> target daily sum (station-like),
    then rescales all intra-day slots proportionally.
    """
    if not daily_models:
        return calib_mm

    for i0, i1 in _build_day_slices(dts_arr):
        raw_day = raw_mm[i0:i1].sum(axis=0)
        calib_day = calib_mm[i0:i1].sum(axis=0)

        for wmo, mask in station_masks.items():
            params = daily_models.get(wmo)
            if params is None:
                continue

            raw_vals = raw_day[mask]
            if raw_vals.size == 0:
                continue

            q_sat, q_st, slope, p_th = params
            target_daily = np.maximum(apply_qm(raw_vals, q_sat, q_st, slope, p_th), 0.0)
            current_daily = calib_day[mask]

            factor = np.ones_like(current_daily, dtype=np.float64)
            pos = current_daily > 0

            # Needed ratio to reach target daily sum.
            needed = np.ones_like(current_daily, dtype=np.float64)
            needed[pos] = target_daily[pos] / current_daily[pos]

            # If current daily sum is zero but target > 0, use raw daily as fallback base.
            miss = (~pos) & (target_daily > 0) & (raw_vals > 0)
            needed[miss] = target_daily[miss] / raw_vals[miss]

            # Apply only where mismatch is meaningful; keep correction soft.
            mismatch = np.abs(needed - 1.0)
            active = mismatch > 0.20
            factor[active] = 1.0 + 0.35 * (needed[active] - 1.0)
            factor = np.clip(factor, 0.85, 1.15)
            calib_mm[i0:i1, mask] *= factor

    return calib_mm


def apply_year_anchor_constraint(calib_mm: np.ndarray,
                                 raw_mm: np.ndarray,
                                 station_masks: dict,
                                 annual_year_ratio: dict,
                                 year: int) -> np.ndarray:
    """
    Year-specific annual anchor from station/satellite annual ratio:
    for each station-zone in a given year nudges annual total to the
    corresponding station-observed annual ratio (if full-year anchor exists).
    """
    if not annual_year_ratio:
        return calib_mm

    raw_annual = raw_mm.sum(axis=0)
    calib_annual = calib_mm.sum(axis=0)

    for wmo, mask in station_masks.items():
        ymap = annual_year_ratio.get(wmo)
        if not ymap:
            continue

        ratio = ymap.get(int(year))
        if ratio is None:
            continue

        raw_vals = raw_annual[mask]
        corr_vals = calib_annual[mask]
        if raw_vals.size == 0:
            continue

        target_vals = np.maximum(raw_vals, 0.0) * float(ratio)

        needed = np.ones_like(corr_vals, dtype=np.float64)
        pos = corr_vals > 0
        needed[pos] = target_vals[pos] / corr_vals[pos]

        miss = (~pos) & (target_vals > 0) & (raw_vals > 0)
        needed[miss] = target_vals[miss] / raw_vals[miss]

        mismatch = np.abs(needed - 1.0)
        active = mismatch > 0.03

        # Strong anchor to enforce interannual realism in extreme years.
        factor = np.ones_like(corr_vals, dtype=np.float64)
        factor[active] = 1.0 + 0.65 * (needed[active] - 1.0)
        factor = np.clip(factor, 0.70, 1.40)
        calib_mm[:, mask] *= factor

    return calib_mm


def _fit_annual_transfer_model(annual_df: pd.DataFrame) -> dict | None:
    """
    Fits station-level annual transfer:
      raw annual total -> station annual total.
    Uses a monotonic empirical mapping over annual quantiles and
    robust linear tails for extrapolation.
    """
    if annual_df is None or annual_df.empty:
        return None

    adf = annual_df.copy()
    adf = adf[
        np.isfinite(adf['P_sat_mm'].values)
        & np.isfinite(adf['P_station_mm'].values)
        & (adf['P_sat_mm'].values > 0)
    ]
    if len(adf) < 6:
        return None

    x = adf['P_sat_mm'].to_numpy(dtype=np.float64)
    y = adf['P_station_mm'].to_numpy(dtype=np.float64)
    ratio_med = float(np.median(np.clip(y / np.maximum(x, 1e-6), 0.2, 5.0)))

    q_levels = np.linspace(0.05, 0.95, 19)
    xq = np.quantile(x, q_levels)
    yq = np.quantile(y, q_levels)

    # Ensure strictly increasing x-grid for interpolation.
    xq_u, idx_u = np.unique(xq, return_index=True)
    yq_u = yq[idx_u]

    if len(xq_u) < 4:
        return {
            'kind': 'ratio',
            'ratio': ratio_med,
            'raw_p10': float(np.percentile(x, 10)),
            'raw_p50': float(np.percentile(x, 50)),
            'raw_p90': float(np.percentile(x, 90)),
        }

    # Keep monotonic response.
    yq_u = np.maximum.accumulate(yq_u)

    n_tail = min(4, len(xq_u))
    dx_low = xq_u[n_tail - 1] - xq_u[0]
    dx_high = xq_u[-1] - xq_u[-n_tail]
    slope_low = (
        (yq_u[n_tail - 1] - yq_u[0]) / dx_low if dx_low > 0 else ratio_med
    )
    slope_high = (
        (yq_u[-1] - yq_u[-n_tail]) / dx_high if dx_high > 0 else ratio_med
    )
    slope_low = float(np.clip(slope_low, 0.1, 6.0))
    slope_high = float(np.clip(slope_high, 0.1, 6.0))

    return {
        'kind': 'empirical',
        'xq': xq_u,
        'yq': yq_u,
        'slope_low': slope_low,
        'slope_high': slope_high,
        'raw_p10': float(np.percentile(x, 10)),
        'raw_p50': float(np.percentile(x, 50)),
        'raw_p90': float(np.percentile(x, 90)),
    }


def _map_annual_raw_to_target(raw_vals: np.ndarray, model: dict) -> np.ndarray:
    """Maps raw annual totals to target annual totals using fitted model."""
    raw = np.asarray(raw_vals, dtype=np.float64)
    out = np.zeros_like(raw, dtype=np.float64)

    if model.get('kind') == 'ratio':
        ratio = float(model.get('ratio', 1.0))
        out = raw * ratio
        out[raw <= 0] = 0.0
        return np.maximum(out, 0.0)

    xq = np.asarray(model['xq'], dtype=np.float64)
    yq = np.asarray(model['yq'], dtype=np.float64)
    slope_low = float(model.get('slope_low', 1.0))
    slope_high = float(model.get('slope_high', 1.0))

    out = np.interp(raw, xq, yq)
    low = raw < xq[0]
    high = raw > xq[-1]
    if np.any(low):
        out[low] = yq[0] + slope_low * (raw[low] - xq[0])
    if np.any(high):
        out[high] = yq[-1] + slope_high * (raw[high] - xq[-1])

    out[raw <= 0] = 0.0
    out[~np.isfinite(out)] = 0.0
    return np.maximum(out, 0.0)


def apply_annual_transfer_constraint(calib_mm: np.ndarray,
                                     raw_mm: np.ndarray,
                                     station_masks: dict,
                                     annual_transfer_models: dict) -> np.ndarray:
    """
    Interannual amplitude control:
    for each station-zone maps raw annual total -> target annual total.
    Correction is adaptive: stronger in extreme dry/wet years.
    """
    if not annual_transfer_models:
        return calib_mm

    raw_annual = raw_mm.sum(axis=0)
    calib_annual = calib_mm.sum(axis=0)

    for wmo, mask in station_masks.items():
        model = annual_transfer_models.get(wmo)
        if model is None:
            continue

        raw_vals = raw_annual[mask]
        corr_vals = calib_annual[mask]
        if raw_vals.size == 0:
            continue

        target_vals = _map_annual_raw_to_target(raw_vals, model)

        needed = np.ones_like(corr_vals, dtype=np.float64)
        pos = corr_vals > 0
        needed[pos] = target_vals[pos] / corr_vals[pos]

        miss = (~pos) & (target_vals > 0) & (raw_vals > 0)
        needed[miss] = target_vals[miss] / raw_vals[miss]

        spread = max(float(model['raw_p90'] - model['raw_p10']), 1.0)
        z = np.abs(raw_vals - float(model['raw_p50'])) / spread
        # Base gain in normal years + stronger gain in extremes.
        gain = 0.25 + 0.35 * np.clip(z, 0.0, 2.0) / 2.0

        mismatch = np.abs(needed - 1.0)
        active = mismatch > 0.08
        factor = np.ones_like(corr_vals, dtype=np.float64)
        factor[active] = 1.0 + gain[active] * (needed[active] - 1.0)

        # Keep physically plausible annual adjustment.
        factor = np.clip(factor, 0.75, 1.30)
        calib_mm[:, mask] *= factor

    return calib_mm


def apply_annual_sanity_guard(calib_mm: np.ndarray,
                              raw_mm: np.ndarray,
                              station_masks: dict,
                              annual_stats: dict) -> np.ndarray:
    """
    Annual sanity guard: keeps calibrated annual totals within station-trained
    annual ratio envelopes (relative to raw annual amount).
    """
    if not annual_stats:
        return calib_mm

    raw_annual = raw_mm.sum(axis=0)
    calib_annual = calib_mm.sum(axis=0)

    for wmo, mask in station_masks.items():
        stats = annual_stats.get(wmo)
        if stats is None:
            continue

        raw_vals = raw_annual[mask]
        corr_vals = calib_annual[mask]
        if raw_vals.size == 0:
            continue

        low_ratio = max(0.0, stats['ratio_p10'] * 0.9)
        high_ratio = max(low_ratio + 1e-6, stats['ratio_p90'] * 1.1)

        lower = raw_vals * low_ratio
        upper = raw_vals * high_ratio

        abs_upper = stats['station_p90'] * 1.2
        if np.isfinite(abs_upper):
            upper = np.minimum(upper, abs_upper)

        factor = np.ones_like(corr_vals, dtype=np.float64)
        over = corr_vals > upper
        needed_over = np.ones_like(corr_vals, dtype=np.float64)
        needed_over[over] = upper[over] / corr_vals[over]

        under = (corr_vals < lower) & (corr_vals > 0)
        needed_under = np.ones_like(corr_vals, dtype=np.float64)
        needed_under[under] = lower[under] / corr_vals[under]

        needed = np.ones_like(corr_vals, dtype=np.float64)
        needed[over] = needed_over[over]
        needed[under] = needed_under[under]

        mismatch = np.abs(needed - 1.0)
        active = mismatch > 0.15
        factor[active] = 1.0 + 0.45 * (needed[active] - 1.0)

        factor = np.clip(factor, 0.85, 1.15)
        calib_mm[:, mask] *= factor

    return calib_mm


# ==================================================================
# 5. Основная функция обработки одного года
# ==================================================================

def process_year(year: int,
                 tif_paths: list,
                 out_dir: str,
                 out_fname: str,
                 nearest_wmo_grid: np.ndarray,
                 seasonal_models: dict,
                 blend_alpha: dict,
                 daily_models: dict,
                 annual_year_ratio: dict,
                 annual_transfer_models: dict,
                 annual_stats: dict,
                 step_min: int,
                 slots_per_3h: int) -> None:
    """
    Merges quarterly rasters for one year, applies:
      3h QM -> disaggregation -> daily volume control
      -> year-anchor annual control -> annual transfer control
      -> annual sanity guard,
    then writes annual calibrated GeoTIFF.

    Parameters
    ---------
    step_min    : source temporal step in minutes (30 for IMERG, 60 for ERA5)
    slots_per_3h: slots per one 3h window (6 for IMERG, 3 for ERA5)
    """
    step_h = step_min / 60.0

    print(f"\n{'='*60}")
    print(f"  Year {year}: {len(tif_paths)} quarter file(s)")
    print(f"{'='*60}")

    all_raw: list = []
    all_dts: list = []
    ref_ds = None

    for tif_path in sorted(tif_paths):
        print(f"  Reading: {os.path.basename(tif_path)}")
        ds = rxr.open_rasterio(tif_path, masked=True)

        data = ds.values.astype(np.float64)  # (bands, rows, cols), mm/h
        data = np.where(np.isnan(data), 0.0, data)
        data = np.maximum(data, 0.0)

        long_name = ds.attrs.get('long_name', None)
        if long_name is None:
            raise ValueError(f"No long_name attribute in {os.path.basename(tif_path)}")

        dts = parse_band_datetimes(long_name)
        for i, dt in enumerate(dts):
            if dt is None:
                print(f"    WARN: failed to parse datetime for band {i}")
                continue
            all_raw.append(data[i])
            all_dts.append(dt)

        if ref_ds is None:
            ref_ds = ds

    if not all_raw:
        print(f"  WARN: no data for year {year}, skipping.")
        return

    # Time sort
    order = np.argsort([dt.value for dt in all_dts])
    raw_arr = np.stack([all_raw[i] for i in order], axis=0)  # (N, R, C), mm/h
    dts_arr = [all_dts[i] for i in order]

    n_slots = raw_arr.shape[0]
    n_rows, n_cols = raw_arr.shape[1], raw_arr.shape[2]

    print(f"  Total slots: {n_slots} ({dts_arr[0]} .. {dts_arr[-1]}), step {step_min} min")

    station_masks = {
        int(wmo): (nearest_wmo_grid == wmo)
        for wmo in np.unique(nearest_wmo_grid)
    }

    # ------------------------------------------------------------------
    # 3h aggregation -> QM -> disaggregation
    # ------------------------------------------------------------------
    calib_arr = np.copy(raw_arr)  # output in mm/h

    n_windows = n_slots // slots_per_3h
    remainder = n_slots % slots_per_3h

    for w in tqdm(range(n_windows), desc=f"  Calibration {year} (3h windows)", leave=False):
        i0 = w * slots_per_3h
        i1 = i0 + slots_per_3h

        window_mmh = raw_arr[i0:i1]          # (slots_per_3h, R, C), mm/h
        window_mm = window_mmh * step_h      # accumulated mm per slot
        sum_3h_mm = window_mm.sum(axis=0)    # 3h total (R, C), mm

        mid_idx = i0 + slots_per_3h // 2
        season = get_season(dts_arr[mid_idx].month)

        calib_3h_mm = apply_qm_to_3h_window(
            sum_3h_mm,
            station_masks,
            season,
            seasonal_models,
            blend_alpha=blend_alpha,
        )

        # Proportional disaggregation back to source slots.
        calib_window_mm = np.zeros_like(window_mm)
        mask_gt0 = sum_3h_mm > 0

        for j in range(slots_per_3h):
            calib_window_mm[j, mask_gt0] = (
                calib_3h_mm[mask_gt0]
                * (window_mm[j, mask_gt0] / sum_3h_mm[mask_gt0])
            )

        calib_arr[i0:i1] = calib_window_mm / step_h

    if remainder > 0:
        print(f"  INFO: {remainder} tail slot(s) left uncorrected (incomplete 3h window).")

    # ------------------------------------------------------------------
    # Daily + Annual guards in accumulated units (mm)
    # ------------------------------------------------------------------
    raw_mm = raw_arr * step_h
    calib_mm = np.maximum(calib_arr, 0.0) * step_h

    if daily_models:
        calib_mm = apply_daily_qm_constraint(
            calib_mm=calib_mm,
            raw_mm=raw_mm,
            dts_arr=dts_arr,
            station_masks=station_masks,
            daily_models=daily_models,
        )

    if annual_year_ratio:
        calib_mm = apply_year_anchor_constraint(
            calib_mm=calib_mm,
            raw_mm=raw_mm,
            station_masks=station_masks,
            annual_year_ratio=annual_year_ratio,
            year=year,
        )

    if annual_transfer_models:
        calib_mm = apply_annual_transfer_constraint(
            calib_mm=calib_mm,
            raw_mm=raw_mm,
            station_masks=station_masks,
            annual_transfer_models=annual_transfer_models,
        )

    if annual_stats:
        calib_mm = apply_annual_sanity_guard(
            calib_mm=calib_mm,
            raw_mm=raw_mm,
            station_masks=station_masks,
            annual_stats=annual_stats,
        )

    calib_arr = np.maximum(calib_mm / step_h, 0.0).astype(np.float32)

    # ------------------------------------------------------------------
    # Write GeoTIFF
    # ------------------------------------------------------------------
    out_path = os.path.join(out_dir, out_fname)
    transform = ref_ds.rio.transform()
    crs = ref_ds.rio.crs

    band_names = str(tuple(
        f"{i}_P_{dt.strftime('%Y%m%d_%H%M')}" for i, dt in enumerate(dts_arr)
    ))

    with rasterio.open(
        out_path, 'w',
        driver='GTiff',
        height=n_rows,
        width=n_cols,
        count=n_slots,
        dtype='float32',
        crs=crs,
        transform=transform,
        compress='lzw',
        predictor=2,
        tiled=True,
        blockxsize=256,
        blockysize=256,
    ) as dst:
        for b in range(n_slots):
            dst.write(calib_arr[b], b + 1)
        dst.update_tags(long_name=band_names)

    size_mb = os.path.getsize(out_path) / 1024 / 1024
    print(f"  Saved: {out_fname} ({size_mb:.1f} MB)")


# ==================================================================
# 6. Главная точка входа
# ==================================================================

def main():
    BASE_DIR = r'd:\Cache\Yandex.Disk\РАЗРАБОТКА\code\imerg2meteo_calib'

    parser = argparse.ArgumentParser(
        description='Применение QM-калибровки к годовым растрам IMERG или ERA5-Land.'
    )
    parser.add_argument(
        '--dataset', required=True, choices=['imerg', 'era5land'],
        help='Набор данных: imerg или era5land'
    )
    parser.add_argument('--zip',   default=None,
                        help='Путь к ZIP-архиву (по умолчанию — последний архив в РНФ25-28/Осадки)')
    parser.add_argument('--out',   default=None,
                        help='Каталог для записи GeoTIFF (по умолчанию из конфига датасета)')
    parser.add_argument('--calib', default=None,
                        help='Каталог с *_calib.csv (по умолчанию из конфига датасета)')
    parser.add_argument('--meteo',
                        default=os.path.join(BASE_DIR, 'data', 'meteo', 'срочные данные_осадки'),
                        help='Каталог с CSV метеостанций')
    parser.add_argument('--year', type=int, default=None,
                        help='Обработать только один год')
    parser.add_argument('--train-start', default=None, dest='train_start',
                        help='Начало тренировочного периода QM (ГГГГ-ММ-ДД)')
    parser.add_argument('--train-end',   default=None, dest='train_end',
                        help='Конец тренировочного периода QM (ГГГГ-ММ-ДД)')
    parser.add_argument('--tmpdir', default=None,
                        help='Временный каталог для распаковки (по умолчанию %%TEMP%%)')
    args = parser.parse_args()

    cfg = DATASET_CONFIG[args.dataset]

    # Пути по умолчанию
    ARCHIVE_DIR = r'd:\Cache\Yandex.Disk\РНФ25-28\Осадки'
    if args.zip is None:
        # Ищем последний подходящий архив
        if args.dataset == 'imerg':
            pattern = os.path.join(ARCHIVE_DIR, 'IMERG_RFACTOR_ANNUAL*.zip')
        else:
            pattern = os.path.join(ARCHIVE_DIR, 'ERA5LAND_RFACTOR_ANNUAL*.zip')
        zips = sorted(glob.glob(pattern))
        if not zips:
            parser.error(f'Не найден архив по шаблону {pattern}. '
                         f'Укажите явно через --zip.')
        args.zip = zips[-1]

    calib_dir = args.calib or os.path.join(BASE_DIR, cfg['calib_subdir'])
    out_dir   = args.out   or os.path.join(BASE_DIR, cfg['out_subdir'])
    os.makedirs(out_dir, exist_ok=True)

    train_start = args.train_start or cfg['train_start']
    train_end   = args.train_end   or cfg['train_end']

    # ------------------------------------------------------------------
    # A) Распаковка архива
    # ------------------------------------------------------------------
    print(f'\n[{args.dataset.upper()}] Архив: {os.path.basename(args.zip)}')
    print(f'[1/4] Распаковка...')
    tmp_root = tempfile.mkdtemp(prefix=f'{args.dataset}_qm_', dir=args.tmpdir)

    try:
        with zipfile.ZipFile(args.zip, 'r') as z:
            z.extractall(tmp_root)
        print(f'  Распаковано в: {tmp_root}')

        # ------------------------------------------------------------------
        # B) Сбор TIF и группировка по годам
        # ------------------------------------------------------------------
        all_tifs = glob.glob(os.path.join(tmp_root, '**', '*.tif'), recursive=True)
        print(f'\n[2/4] Найдено TIF-файлов: {len(all_tifs)}')

        year_tifs = defaultdict(list)
        for tif in all_tifs:
            m = cfg['file_re'].search(os.path.basename(tif))
            if m:
                year_tifs[int(m.group(1))].append(tif)
            else:
                print(f'  WARN: не распознан файл {os.path.basename(tif)}')

        if not year_tifs:
            print('ERROR: Ни один файл не соответствует шаблону датасета.')
            print(f'  Ожидаемый шаблон: {cfg["file_re"].pattern}')
            sys.exit(1)

        years_available = sorted(year_tifs)
        print(f'  Годы: {years_available[0]}–{years_available[-1]}  '
              f'({len(years_available)} лет)')

        if args.year is not None:
            if args.year not in year_tifs:
                print(f'ERROR: год {args.year} отсутствует в архиве.')
                sys.exit(1)
            years_to_process = [args.year]
        else:
            years_to_process = years_available

        # ------------------------------------------------------------------
        # C) Метаданные станций + KNN + QM-модели
        # ------------------------------------------------------------------
        print(f'\n[3/4] Построение QM-моделей и KNN-сетки...')

        print('  Загрузка метаданных станций...')
        st_meta = load_station_metadata(args.meteo)
        print(f'  Загружено {len(st_meta)} станций.')

        print(f'  Обучение QM-моделей ({train_start} – {train_end})...')
        qm_bundle = precalculate_qm_models(calib_dir, train_start, train_end)

        # KNN по первому TIF из первого года
        sample_tif = year_tifs[years_to_process[0]][0]
        ds_sample  = rxr.open_rasterio(sample_tif)
        lons, lats = np.meshgrid(ds_sample.x.values, ds_sample.y.values)
        pixel_coords = np.column_stack([lons.ravel(), lats.ravel()])

        st_coords = st_meta[['lon', 'lat']].values
        tree = cKDTree(st_coords)
        _, nn_idx = tree.query(pixel_coords)
        nearest_wmo_grid = st_meta['wmo_index'].values[nn_idx].reshape(lons.shape)

        print(f'  KNN-сетка: {nearest_wmo_grid.shape} пикселей → '
              f'{len(np.unique(nearest_wmo_grid))} уникальных станций')

        # ------------------------------------------------------------------
        # D) Основной цикл по годам
        # ------------------------------------------------------------------
        print(f'\n[4/4] Калибровка {len(years_to_process)} год(а/ов)...')

        step_min    = cfg['step_min']
        slots_per_3h = cfg['slots_per_3h']

        for yr in years_to_process:
            out_fname = cfg['out_fname_tpl'].format(year=yr)
            out_path  = os.path.join(out_dir, out_fname)

            if os.path.exists(out_path):
                print(f'\n  → {out_fname} уже существует, пропуск.')
                continue

            process_year(
                year=yr,
                tif_paths=year_tifs[yr],
                out_dir=out_dir,
                out_fname=out_fname,
                nearest_wmo_grid=nearest_wmo_grid,
                seasonal_models=qm_bundle['seasonal'],
                blend_alpha=qm_bundle.get('blend_alpha', {}),
                daily_models=qm_bundle['daily'],
                annual_year_ratio=qm_bundle.get('annual_year_ratio', {}),
                annual_transfer_models=qm_bundle.get('annual_transfer', {}),
                annual_stats=qm_bundle['annual'],
                step_min=step_min,
                slots_per_3h=slots_per_3h,
            )

        print(f'\nГотово. Результаты: {out_dir}')

    finally:
        shutil.rmtree(tmp_root, ignore_errors=True)
        print(f'Временный каталог удалён: {tmp_root}')


if __name__ == '__main__':
    main()
