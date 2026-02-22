"""
apply_qm_to_rasters.py
======================
Применяет разработанную методику QM-калибровки к архиву квартальных
IMERG V07 растров (1 ч, мм/ч), группирует их по годам и сохраняет
годовые откалиброванные GeoTIFF.

Входной архив: IMERG_RFACTOR_ANNUAL-*.zip
  Содержит TIF-файлы вида:
    IMERG_RFACTOR_ANNUAL/IMERG_V07_P1h_mm_YYYY_QN_permanent.tif
  Каждый TIF: (N_hours, n_rows, n_cols), бэнды имеют атрибут long_name:
    tuple of 'IDX_P_YYYYMMDD_HHmm'

Методика:
  - Пространственный KNN: каждый пиксель → ближайшая метеостанция
  - QM-калибровка: Full-Distribution QM (1000 квантилей) + Volume Scaling
    по сезонам (DJF, MAM, JJA, SON) из предварительно обученных моделей
  - Группировка квартальных файлов → годовые выходные файлы

Использование:
  python apply_qm_to_rasters.py [--zip PATH_TO_ZIP] [--out OUTPUT_DIR]
                                 [--calib CALIB_DIR] [--meteo METEO_DIR]
                                 [--year YEAR]       # обработать только один год
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
    from rasterio.transform import from_bounds
except ImportError:
    print("ERROR: Требуются пакеты rioxarray и rasterio.")
    print("  pip install rioxarray rasterio")
    sys.exit(1)

# ------------------------------------------------------------------
# Добавляем src в sys.path чтобы импортировать qm_calibration
# ------------------------------------------------------------------
_SRC = os.path.dirname(os.path.abspath(__file__))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from qm_calibration import get_season, apply_qm


# ==================================================================
# 1.  Загрузка метаданных станций
# ==================================================================

def load_station_metadata(meteo_dir: str) -> pd.DataFrame:
    """
    Читает WMO-индекс, координаты (X=lon, Y=lat) из заголовков CSV.
    """
    csv_files = glob.glob(os.path.join(meteo_dir, "*.csv"))
    meta = []
    for f in csv_files:
        try:
            df = pd.read_csv(f, sep=';', encoding='cp866', nrows=1)
            meta.append({
                'wmo_index': int(df['Index'].iloc[0]),
                'name':      str(df['StationName'].iloc[0]),
                'lon':       float(df['X'].iloc[0]),
                'lat':       float(df['Y'].iloc[0]),
            })
        except Exception:
            pass
    if not meta:
        raise RuntimeError(f"Не найдено ни одной метеостанции в {meteo_dir}")
    return pd.DataFrame(meta)


# ==================================================================
# 2.  Предварительный расчёт QM-моделей из calib CSV
# ==================================================================

def precalculate_qm_models(calib_dir: str,
                            train_start: str = '2001-01-01',
                            train_end:   str = '2015-12-31') -> dict:
    """
    Возвращает dict[season][wmo_index] = (q_sat, q_station, slope_extrap, p_th).

    Модели обучаются на тренировочном периоде (train_start…train_end).
    Используется Full-Distribution QM (1000 квантилей, весь CDF).
    """
    models = {s: {} for s in ('DJF', 'MAM', 'JJA', 'SON')}
    calib_files = glob.glob(os.path.join(calib_dir, "*_calib.csv"))

    if not calib_files:
        raise RuntimeError(f"Не найдено файлов *_calib.csv в {calib_dir}")

    num_q = 1000
    q_levels = np.linspace(0.0, 1.0, num_q + 2)[1:-1]  # без 0 и 1

    for f in tqdm(calib_files, desc="Обучение QM-моделей"):
        try:
            df = pd.read_csv(f)
        except Exception as e:
            print(f"  Пропуск {f}: {e}")
            continue

        if 'wmo_index' not in df.columns:
            continue

        wmo = int(df['wmo_index'].iloc[0])
        df['datetime'] = pd.to_datetime(df['datetime'])
        if 'season' not in df.columns:
            df['season'] = df['datetime'].dt.month.apply(get_season)

        train_mask = (df['datetime'] >= train_start) & (df['datetime'] <= train_end)
        train_df   = df[train_mask]

        for season in models:
            sd = train_df[train_df['season'] == season]
            p_sat = sd['P_sat_mm'].values
            p_st  = sd['P_station_mm'].values

            if len(p_sat) < 30 or len(p_st) < 30:
                continue

            q_sat = np.quantile(p_sat, q_levels)
            q_st  = np.quantile(p_st,  q_levels)

            # Экстраполяция правого хвоста
            if (q_sat[-1] - q_sat[-5]) > 0:
                slope = (q_st[-1] - q_st[-5]) / (q_sat[-1] - q_sat[-5])
            else:
                slope = 1.0

            # p_th = 0 (Full-Distribution — весь CDF, не только wet-days)
            models[season][wmo] = (q_sat, q_st, slope, 0.0)

    total = sum(len(v) for v in models.values())
    print(f"  Готово. Всего пар сезон×станция: {total}")
    return models


# ==================================================================
# 3.  Разбор имён бэндов → datetime
# ==================================================================

def parse_band_datetimes(long_name_tuple) -> list:
    """
    Разбирает кортеж long_name вида ('0_P_20010101_0000', ...) → list[datetime].
    """
    dts = []
    pat = re.compile(r'P_(\d{8})_(\d{4})')
    for name in long_name_tuple:
        m = pat.search(str(name))
        if m:
            dt = pd.to_datetime(m.group(1) + m.group(2), format='%Y%m%d%H%M')
            dts.append(dt)
        else:
            dts.append(None)
    return dts


# ==================================================================
# 4.  Применение QM к одному массиву пикселей одного часа
# ==================================================================

def apply_qm_to_hour(raw_hour: np.ndarray,
                     nearest_wmo_grid: np.ndarray,
                     season: str,
                     models: dict) -> np.ndarray:
    """
    Применяет QM-коррекцию к 2-D массиву одного часа.
    Пиксели без модели для соответствующей станции остаются без изменений.

    raw_hour : (rows, cols), мм/ч
    Возвращает откалиброванный массив той же формы.
    """
    calib = np.copy(raw_hour).astype(float)
    season_models = models.get(season, {})

    unique_wmos = np.unique(nearest_wmo_grid)
    for wmo in unique_wmos:
        pixel_mask = nearest_wmo_grid == wmo
        vals = raw_hour[pixel_mask]

        if len(vals) == 0:
            continue

        params = season_models.get(wmo, None)
        if params is None:
            continue  # нет модели — оставляем как есть

        q_sat, q_st, slope, p_th = params
        corrected = apply_qm(vals, q_sat, q_st, slope, p_th)
        corrected  = np.maximum(corrected, 0.0)
        calib[pixel_mask] = corrected

    return calib


# ==================================================================
# 5.  Основная функция обработки одного года
# ==================================================================

def process_year(year: int,
                 tif_paths: list,
                 out_dir: str,
                 nearest_wmo_grid: np.ndarray,
                 models: dict,
                 reference_ds=None) -> None:
    """
    Обрабатывает все квартальные TIF одного года, объединяет их
    и сохраняет годовой откалиброванный GeoTIFF.

    tif_paths : список путей к TIF-файлам за данный год (Q1..Q4), 
                уже упорядоченных по дате.
    """
    print(f"\n{'='*60}")
    print(f"  Год {year}: {len(tif_paths)} квартал(а/ов)")
    print(f"{'='*60}")

    all_bands_raw  = []   # list of 2D arrays (одночасовые)
    all_datetimes  = []   # list of pd.Timestamp

    for tif_path in sorted(tif_paths):
        fname = os.path.basename(tif_path)
        print(f"  Чтение: {fname}")
        ds = rxr.open_rasterio(tif_path, masked=True)

        # Заменяем NaN → 0 (IMERG не имеет NoData в данном регионе)
        data = ds.values  # (bands, rows, cols)
        data = np.where(np.isnan(data), 0.0, data)
        data = np.maximum(data, 0.0)

        # Разбираем дата-время бэндов
        long_name = ds.attrs.get('long_name', None)
        if long_name is None:
            raise ValueError(f"Нет атрибута long_name в {fname}")

        datetimes = parse_band_datetimes(long_name)

        for i, dt in enumerate(datetimes):
            if dt is None:
                print(f"    WARN: не удалось разобрать дату бэнда {i}")
                continue
            all_bands_raw.append(data[i, :, :])
            all_datetimes.append(dt)

        # Сохраняем референсный датасет для геопривязки
        if reference_ds is None:
            reference_ds = ds

    if not all_bands_raw:
        print(f"  WARN: Нет данных для года {year}, пропуск.")
        return

    # Сортировка по времени
    order = np.argsort([dt.value for dt in all_datetimes])
    all_bands_sorted = [all_bands_raw[i] for i in order]
    dts_sorted       = [all_datetimes[i] for i in order]

    n_hours = len(all_bands_sorted)
    print(f"  Всего часов: {n_hours}  ({dts_sorted[0]} … {dts_sorted[-1]})")

    # Калибровка почасовых бэндов
    calibrated_bands = []
    for i, (raw_h, dt) in enumerate(
            tqdm(zip(all_bands_sorted, dts_sorted),
                 total=n_hours, desc=f"  Калибровка {year}", leave=False)):

        season = get_season(dt.month)
        calib_h = apply_qm_to_hour(raw_h, nearest_wmo_grid, season, models)
        calibrated_bands.append(calib_h)

    # Сборка и запись GeoTIFF
    out_array = np.stack(calibrated_bands, axis=0).astype(np.float32)
    # out_array shape: (n_hours, rows, cols)

    out_fname = f"IMERG_V07_P1h_mm_{year}_calib_qm.tif"
    out_path  = os.path.join(out_dir, out_fname)

    # Строим long_name для выходного файла
    band_names = tuple(
        f"{i}_P_{dt.strftime('%Y%m%d_%H%M')}" for i, dt in enumerate(dts_sorted)
    )

    # Сохраняем через rasterio напрямую
    ref = reference_ds
    transform = ref.rio.transform()
    crs        = ref.rio.crs

    with rasterio.open(
        out_path, 'w',
        driver  = 'GTiff',
        height  = out_array.shape[1],
        width   = out_array.shape[2],
        count   = out_array.shape[0],
        dtype   = 'float32',
        crs     = crs,
        transform=transform,
        compress='lzw',
        predictor=2,
        tiled   = True,
        blockxsize=256,
        blockysize=256,
    ) as dst:
        for b in range(out_array.shape[0]):
            dst.write(out_array[b], b + 1)
        # Записываем имена бэндов через descriptions
        dst.update_tags(long_name=str(band_names))

    size_mb = os.path.getsize(out_path) / 1024 / 1024
    print(f"  ✓ Сохранено: {out_fname}  ({size_mb:.1f} МБ)")


# ==================================================================
# 6.  Главная точка входа
# ==================================================================

def main():
    BASE_DIR  = r"d:\Cache\Yandex.Disk\РАЗРАБОТКА\code\imerg2meteo_calib"
    ZIP_DEFAULT = r"d:\Cache\Yandex.Disk\РНФ25-28\Осадки\IMERG_RFACTOR_ANNUAL-20260222T084003Z-1-001.zip"

    parser = argparse.ArgumentParser(
        description="Применение QM-калибровки к годовым IMERG растрам из архива."
    )
    parser.add_argument('--zip',   default=ZIP_DEFAULT,
                        help="Путь к ZIP-архиву с квартальными TIF")
    parser.add_argument('--out',   default=os.path.join(BASE_DIR, 'output', 'imerg_rfactor_calib'),
                        help="Каталог для записи откалиброванных GeoTIFF")
    parser.add_argument('--calib', default=os.path.join(BASE_DIR, 'output', 'calib_imerg'),
                        help="Каталог с файлами *_calib.csv (результаты main.py --dataset imerg)")
    parser.add_argument('--meteo', default=os.path.join(BASE_DIR, 'data', 'meteo', 'срочные данные_осадки'),
                        help="Каталог с CSV метеостанций (для координат)")
    parser.add_argument('--year',  type=int, default=None,
                        help="Обработать только один год (например, --year 2010)")
    parser.add_argument('--train-start', default='2001-01-01', dest='train_start',
                        help="Начало тренировочного периода QM (ГГГГ-ММ-ДД)")
    parser.add_argument('--train-end',   default='2015-12-31', dest='train_end',
                        help="Конец тренировочного периода QM (ГГГГ-ММ-ДД)")
    parser.add_argument('--tmpdir', default=None,
                        help="Временный каталог для распаковки TIF (по умолчанию %TEMP%)")
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)

    # ------------------------------------------------------------------
    # A) Распаковка архива во временный каталог
    # ------------------------------------------------------------------
    print(f"\n[1/4] Распаковка архива: {os.path.basename(args.zip)}")
    tmp_root = tempfile.mkdtemp(prefix='imerg_qm_', dir=args.tmpdir)
    try:
        with zipfile.ZipFile(args.zip, 'r') as z:
            z.extractall(tmp_root)
        print(f"  Распаковано в: {tmp_root}")

        # ------------------------------------------------------------------
        # B) Сбор TIF-файлов и группировка по годам
        # ------------------------------------------------------------------
        tif_pattern = os.path.join(tmp_root, '**', '*.tif')
        all_tifs = glob.glob(tif_pattern, recursive=True)
        print(f"\n[2/4] Найдено TIF-файлов: {len(all_tifs)}")

        # Шаблон: IMERG_V07_P1h_mm_YYYY_QN_permanent.tif
        tif_re = re.compile(r'IMERG_V07_P1h_mm_(\d{4})_Q(\d)_permanent\.tif', re.IGNORECASE)

        year_tifs = defaultdict(list)
        for tif in all_tifs:
            m = tif_re.search(os.path.basename(tif))
            if m:
                yr = int(m.group(1))
                year_tifs[yr].append(tif)
            else:
                print(f"  WARN: не распознан файл {os.path.basename(tif)}")

        years_available = sorted(year_tifs.keys())
        print(f"  Годы в архиве: {years_available[0]}–{years_available[-1]}  "
              f"(всего {len(years_available)} лет)")

        if args.year is not None:
            if args.year not in year_tifs:
                print(f"ERROR: год {args.year} отсутствует в архиве.")
                sys.exit(1)
            years_to_process = [args.year]
        else:
            years_to_process = years_available

        # ------------------------------------------------------------------
        # C) Метаданные станций + KNN-сетка (строится один раз)
        # ------------------------------------------------------------------
        print(f"\n[3/4] Построение QM-моделей и KNN-сетки...")

        print("  Загрузка метаданных станций...")
        st_meta = load_station_metadata(args.meteo)
        print(f"  Загружено {len(st_meta)} станций.")

        print("  Предварительный расчёт QM-моделей...")
        qm_models = precalculate_qm_models(
            args.calib,
            train_start=args.train_start,
            train_end=args.train_end,
        )

        # Определяем пространственную сетку растров из первого TIF
        sample_tif = year_tifs[years_to_process[0]][0]
        ds_sample  = rxr.open_rasterio(sample_tif)
        lons_grid, lats_grid = np.meshgrid(ds_sample.x.values, ds_sample.y.values)
        pixel_coords = np.column_stack([lons_grid.ravel(), lats_grid.ravel()])

        st_coords = st_meta[['lon', 'lat']].values
        tree = cKDTree(st_coords)
        _, nearest_idx = tree.query(pixel_coords)
        nearest_wmo_flat = st_meta['wmo_index'].values[nearest_idx]
        nearest_wmo_grid = nearest_wmo_flat.reshape(lons_grid.shape)

        print(f"  KNN-сетка: {nearest_wmo_grid.shape} пикселей → "
              f"{len(np.unique(nearest_wmo_grid))} уникальных станций")

        # ------------------------------------------------------------------
        # D) Основной цикл по годам
        # ------------------------------------------------------------------
        print(f"\n[4/4] Калибровка {len(years_to_process)} год(а/ов)...")

        for yr in years_to_process:
            # Пропускаем уже обработанные
            out_fname = f"IMERG_V07_P1h_mm_{yr}_calib_qm.tif"
            out_path  = os.path.join(args.out, out_fname)
            if os.path.exists(out_path):
                print(f"\n  → {out_fname} уже существует, пропуск.")
                continue

            process_year(
                year=yr,
                tif_paths=year_tifs[yr],
                out_dir=args.out,
                nearest_wmo_grid=nearest_wmo_grid,
                models=qm_models,
            )

        print(f"\n✓ Готово. Результаты сохранены в: {args.out}")

    finally:
        # Удаляем временный каталог
        shutil.rmtree(tmp_root, ignore_errors=True)
        print(f"  Временный каталог удалён: {tmp_root}")


if __name__ == '__main__':
    main()
