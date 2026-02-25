"""
Compare raw IMERG raster pixel at Kazan vs point-extracted P_sat_mm from calib CSV.
"""
import sys
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
import zipfile, tempfile, os, re, shutil
import numpy as np, pandas as pd
import rasterio  # используем rasterio напрямую — cleanly закрывается

ZIP = r'd:\Cache\Yandex.Disk\РНФ25-28\Осадки\IMERG_RFACTOR_ANNUAL-20260223T120530Z-1-001.zip'
CALIB_CSV = r'd:\Cache\Yandex.Disk\РАЗРАБОТКА\code\imerg2meteo_calib\output\calib_imerg\Казань_27595_calib.csv'
KAZAN_LON, KAZAN_LAT = 49.2, 55.73333333
YEAR = 2001

pat = re.compile(r'P_(\d{8})_(\d{4})')

# 1. Распакуем нужные кварталы во временную папку
z = zipfile.ZipFile(ZIP)
quarters = sorted([n for n in z.namelist() if f'_{YEAR}_Q' in n])
print(f"Кварталов за {YEAR}: {len(quarters)}")

tmpdir = tempfile.mkdtemp()
try:
    for qfile in quarters:
        tmp = os.path.join(tmpdir, os.path.basename(qfile))
        with open(tmp, 'wb') as f:
            f.write(z.read(qfile))

    # 2. Определяем пиксель Казани по первому файлу
    first_tif = os.path.join(tmpdir, os.path.basename(quarters[0]))
    with rasterio.open(first_tif) as src:
        # Преобразуем координаты в индексы пикселя
        row, col = src.index(KAZAN_LON, KAZAN_LAT)
        lon_center = src.xy(row, col)[0]
        lat_center = src.xy(row, col)[1]
        print(f"Kazan pixel: row={row} col={col}  "
              f"lon={lon_center:.4f} lat={lat_center:.4f}")
        n_bands_q1 = src.count
        tags = src.tags()

    # 3. Читаем временны́е ряды из всех кварталов
    raw_slots_mmh = []
    raw_dts = []

    for qfile in quarters:
        tmp = os.path.join(tmpdir, os.path.basename(qfile))
        with rasterio.open(tmp) as src:
            tags = src.tags()
            long_name_str = tags.get('long_name', '')
            # Парсим tuple из строки
            try:
                band_names = eval(long_name_str)
            except Exception:
                band_names = [str(i) for i in range(src.count)]

            # Читаем ОДИН пиксель для всех бэндов
            # rasterio.sample возвращает значения по (col, row) — (x, y)
            data_1d = np.array([
                src.read(b + 1)[row, col] for b in range(src.count)
            ], dtype=np.float64)

        for i, name in enumerate(band_names):
            m = pat.search(str(name))
            if m:
                dt = pd.to_datetime(m.group(1) + m.group(2), format='%Y%m%d%H%M')
                val = max(float(data_1d[i]), 0.0)
                raw_slots_mmh.append(val)
                raw_dts.append(dt)

finally:
    shutil.rmtree(tmpdir, ignore_errors=True)

raw_slots_mmh = np.array(raw_slots_mmh)  # mm/h intensity, 30-min slots
raw_mm_per_slot = raw_slots_mmh * 0.5    # mm per 30-min step

print(f"\nСырой растр {YEAR}: {len(raw_slots_mmh)} слотов, "
      f"ненулевых={np.sum(raw_slots_mmh > 0)}")
print(f"  Годовая сумма (× 0.5h): {raw_mm_per_slot.sum():.1f} мм")
print(f"  Max intensity: {raw_slots_mmh.max():.3f} мм/ч")

# 3h окна (6 слотов × 30 мин = 3ч)
n3h = len(raw_mm_per_slot) // 6
sums_3h_raw = np.array([raw_mm_per_slot[i*6:(i+1)*6].sum() for i in range(n3h)])
print(f"  3h окон: {n3h}, wet(>0): {np.sum(sums_3h_raw > 0)}")
print(f"  Годовая сумма (3h): {sums_3h_raw.sum():.1f} мм — д.б. равно строке выше")

# 4. Точечные данные из calib CSV
df = pd.read_csv(CALIB_CSV)
df['datetime'] = pd.to_datetime(df['datetime'])
df['year'] = df['datetime'].dt.year
df2001 = df[df['year'] == YEAR]
ps = df2001['P_sat_mm'].values
pst = df2001['P_station_mm'].values

print(f"\nТочечный CSV {YEAR}: {len(df2001)} строк (3h)")
print(f"  P_sat_mm  сумма: {ps.sum():.1f} мм")
print(f"  P_station сумма: {pst.sum():.1f} мм")

# 5. Сравнение квантилей wet-событий
ps_wet  = ps[ps > 0]
raw_wet = sums_3h_raw[sums_3h_raw > 0]

print(f"\nWet 3h events: растр={len(raw_wet)}, точечный={len(ps_wet)}")
print(f"\nКвантили wet 3h (мм/3ч):")
print(f"{'P':>4}  {'Raw растр':>12}  {'Точечный CSV':>14}")
for q in [25, 50, 75, 90, 95, 99]:
    rv = np.percentile(raw_wet, q) if len(raw_wet) else np.nan
    pv = np.percentile(ps_wet, q) if len(ps_wet) else np.nan
    ratio = rv/pv if pv > 0 else np.nan
    print(f"P{q:>2}  {rv:>12.4f}  {pv:>14.4f}  ratio={ratio:.2f}")

print(f"\nВывод:")
print(f"  Raw растр / P_sat_mm (суммы): {raw_mm_per_slot.sum() / ps.sum():.2f}x")
