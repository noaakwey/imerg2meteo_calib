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
        'out_subdir': os.path.join('output', 'imerg_rfactor_calib'),
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

def load_station_metadata(meteo_dir: str) -> pd.DataFrame:
    """Читает WMO-индекс и координаты (X=lon, Y=lat) из заголовков CSV."""
    csv_files = glob.glob(os.path.join(meteo_dir, '*.csv'))
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
        raise RuntimeError(f'Не найдено ни одной метеостанции в {meteo_dir}')
    return pd.DataFrame(meta)


# ==================================================================
# 2. Предварительный расчёт QM-моделей из calib CSV
# ==================================================================

def precalculate_qm_models(calib_dir: str,
                            train_start: str,
                            train_end: str) -> dict:
    """
    Строит Full-Distribution QM (1000 квантилей) сезонные модели.
    Обучение на 3-часовых данных из *_calib.csv (колонки P_sat_mm / P_station_mm).

    Возвращает dict[season][wmo_index] = (q_sat, q_station, slope, p_th=0.0).
    """
    models = {s: {} for s in ('DJF', 'MAM', 'JJA', 'SON')}
    calib_files = glob.glob(os.path.join(calib_dir, '*_calib.csv'))

    if not calib_files:
        raise RuntimeError(f'Не найдено файлов *_calib.csv в {calib_dir}')

    num_q   = 1000
    q_levels = np.linspace(0.0, 1.0, num_q + 2)[1:-1]

    for f in tqdm(calib_files, desc='Обучение QM-моделей'):
        try:
            df = pd.read_csv(f)
        except Exception as e:
            print(f'  Пропуск {f}: {e}')
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
            sd    = train_df[train_df['season'] == season]
            p_sat = sd['P_sat_mm'].values
            p_st  = sd['P_station_mm'].values

            if len(p_sat) < 30 or len(p_st) < 30:
                continue

            q_sat = np.quantile(p_sat, q_levels)
            q_st  = np.quantile(p_st,  q_levels)

            slope = ((q_st[-1] - q_st[-5]) / (q_sat[-1] - q_sat[-5])
                     if (q_sat[-1] - q_sat[-5]) > 0 else 1.0)

            models[season][wmo] = (q_sat, q_st, slope, 0.0)

    total = sum(len(v) for v in models.values())
    print(f'  Готово. Всего пар сезон x станция: {total}')
    return models


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
                           nearest_wmo_grid: np.ndarray,
                           season: str,
                           models: dict) -> np.ndarray:
    """
    Применяет QM-коррекцию к 2D-массиву 3-часовой суммы осадков (мм).
    Пиксели без модели остаются без изменений.
    """
    calib = np.copy(raw_3h_mm).astype(np.float64)
    season_models = models.get(season, {})

    for wmo in np.unique(nearest_wmo_grid):
        mask = nearest_wmo_grid == wmo
        vals = raw_3h_mm[mask]

        if len(vals) == 0:
            continue

        params = season_models.get(wmo, None)
        if params is None:
            continue

        q_sat, q_st, slope, p_th = params
        corrected = apply_qm(vals, q_sat, q_st, slope, p_th)
        corrected  = np.maximum(corrected, 0.0)
        calib[mask] = corrected

    return calib


# ==================================================================
# 5. Основная функция обработки одного года
# ==================================================================

def process_year(year: int,
                 tif_paths: list,
                 out_dir: str,
                 out_fname: str,
                 nearest_wmo_grid: np.ndarray,
                 models: dict,
                 step_min: int,
                 slots_per_3h: int) -> None:
    """
    Объединяет квартальные TIF одного года, применяет
    3h-агрегацию → QM → дезагрегацию и сохраняет годовой GeoTIFF.

    Параметры
    ---------
    step_min    : шаг данных в минутах (30 для IMERG, 60 для ERA5)
    slots_per_3h: слотов в одном 3ч окне (6 или 3 соответственно)
    """
    step_h = step_min / 60.0          # шаг в часах

    print(f"\n{'='*60}")
    print(f"  Год {year}: {len(tif_paths)} квартал(а/ов)")
    print(f"{'='*60}")

    all_raw: list = []     # (rows, cols) float64
    all_dts: list  = []    # pd.Timestamp

    ref_ds = None

    for tif_path in sorted(tif_paths):
        print(f'  Чтение: {os.path.basename(tif_path)}')
        ds = rxr.open_rasterio(tif_path, masked=True)

        data = ds.values.astype(np.float64)    # (bands, rows, cols), мм/ч
        data = np.where(np.isnan(data), 0.0, data)
        data = np.maximum(data, 0.0)

        long_name = ds.attrs.get('long_name', None)
        if long_name is None:
            raise ValueError(f'Нет атрибута long_name в {os.path.basename(tif_path)}')

        dts = parse_band_datetimes(long_name)

        for i, dt in enumerate(dts):
            if dt is None:
                print(f'    WARN: не удалось разобрать дату бэнда {i}')
                continue
            all_raw.append(data[i])
            all_dts.append(dt)

        if ref_ds is None:
            ref_ds = ds

    if not all_raw:
        print(f'  WARN: Нет данных для года {year}, пропуск.')
        return

    # Сортировка по времени
    order     = np.argsort([dt.value for dt in all_dts])
    raw_arr   = np.stack([all_raw[i] for i in order], axis=0)   # (N, R, C) мм/ч
    dts_arr   = [all_dts[i] for i in order]

    n_slots   = raw_arr.shape[0]
    n_rows, n_cols = raw_arr.shape[1], raw_arr.shape[2]

    print(f'  Всего слотов: {n_slots}  '
          f'({dts_arr[0]} … {dts_arr[-1]}), шаг {step_min} мин')

    # ------------------------------------------------------------------
    # 3h-агрегация → QM → дезагрегация
    # ------------------------------------------------------------------
    calib_arr = np.copy(raw_arr)   # выход в тех же единицах мм/ч

    n_windows = n_slots // slots_per_3h
    remainder = n_slots % slots_per_3h

    for w in tqdm(range(n_windows), desc=f'  Калибровка {year} (3h окна)', leave=False):
        i0 = w * slots_per_3h
        i1 = i0 + slots_per_3h

        window_mmh = raw_arr[i0:i1]          # (slots_per_3h, R, C) мм/ч

        # Накопленные мм за каждый слот = интенсивность × шаг_ч
        window_mm  = window_mmh * step_h      # (slots_per_3h, R, C) мм

        # 3ч сумма
        sum_3h_mm  = window_mm.sum(axis=0)    # (R, C) мм

        # Определяем сезон по середине окна
        mid_idx = i0 + slots_per_3h // 2
        season  = get_season(dts_arr[mid_idx].month)

        # QM-коррекция 3ч суммы
        calib_3h_mm = apply_qm_to_3h_window(sum_3h_mm, nearest_wmo_grid,
                                             season, models)   # (R, C) мм

        # Пропорциональная дезагрегация обратно по слотам
        # weight[j] = window_mm[j] / sum_3h_mm  (там где сумма > 0)
        calib_window_mm = np.zeros_like(window_mm)
        mask_gt0 = sum_3h_mm > 0

        for j in range(slots_per_3h):
            calib_window_mm[j, mask_gt0] = (
                calib_3h_mm[mask_gt0]
                * (window_mm[j, mask_gt0] / sum_3h_mm[mask_gt0])
            )
            # Там где 3ч сумма = 0, остаётся 0 (нет осадков)

        # Обратно в мм/ч
        calib_arr[i0:i1] = calib_window_mm / step_h

    # Хвостовые слоты (если N не делится на slots_per_3h) — без коррекции
    if remainder > 0:
        print(f'  INFO: {remainder} хвостовых слот(а) без QM (неполное 3ч окно).')

    calib_arr = np.maximum(calib_arr, 0.0).astype(np.float32)

    # ------------------------------------------------------------------
    # Запись GeoTIFF
    # ------------------------------------------------------------------
    out_path  = os.path.join(out_dir, out_fname)
    transform = ref_ds.rio.transform()
    crs       = ref_ds.rio.crs

    band_names = str(tuple(
        f'{i}_P_{dt.strftime("%Y%m%d_%H%M")}' for i, dt in enumerate(dts_arr)
    ))

    with rasterio.open(
        out_path, 'w',
        driver    = 'GTiff',
        height    = n_rows,
        width     = n_cols,
        count     = n_slots,
        dtype     = 'float32',
        crs       = crs,
        transform = transform,
        compress  = 'lzw',
        predictor = 2,
        tiled     = True,
        blockxsize = 256,
        blockysize = 256,
    ) as dst:
        for b in range(n_slots):
            dst.write(calib_arr[b], b + 1)
        dst.update_tags(long_name=band_names)

    size_mb = os.path.getsize(out_path) / 1024 / 1024
    print(f'  Сохранено: {out_fname}  ({size_mb:.1f} МБ)')


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
        qm_models = precalculate_qm_models(calib_dir, train_start, train_end)

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
                models=qm_models,
                step_min=step_min,
                slots_per_3h=slots_per_3h,
            )

        print(f'\nГотово. Результаты: {out_dir}')

    finally:
        shutil.rmtree(tmp_root, ignore_errors=True)
        print(f'Временный каталог удалён: {tmp_root}')


if __name__ == '__main__':
    main()
