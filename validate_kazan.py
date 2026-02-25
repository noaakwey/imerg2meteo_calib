"""Validation: raster annual sum at Kazan station vs station observations."""
import glob, re
import numpy as np
import pandas as pd
import rioxarray as rxr

# Kazan: X=lon=49.2, Y=lat=55.73 (from meteo CSV)
KAZAN_LON = 49.2
KAZAN_LAT = 55.73333333

RASTER_DIR = r'd:\Cache\Yandex.Disk\Р РђР—Р РђР‘РћРўРљРђ\code\imerg2meteo_calib\output\imerg_rfactor_calib_v5_year_anchor'
CALIB_CSV  = r'd:\Cache\Yandex.Disk\Р РђР—Р РђР‘РћРўРљРђ\code\imerg2meteo_calib\output\calib_imerg\РљР°Р·Р°РЅСЊ_27595_calib.csv'

tifs = sorted(glob.glob(RASTER_DIR + r'\IMERG_V07_P30min_mmh_*_calib_qm.tif'))
print(f"Found {len(tifs)} rasters")

# Index of Kazan pixel in raster grid
ds0 = rxr.open_rasterio(tifs[0])
ix = int(np.argmin(np.abs(ds0.x.values - KAZAN_LON)))
iy = int(np.argmin(np.abs(ds0.y.values - KAZAN_LAT)))
print(f"Kazan pixel: ix={ix} lon={float(ds0.x.values[ix]):.4f}, "
      f"iy={iy} lat={float(ds0.y.values[iy]):.4f}")
print()

rows = []
for tif in tifs:
    yr = int(re.search(r'(\d{4})_calib', tif).group(1))
    ds = rxr.open_rasterio(tif)
    vals = ds.values[:, iy, ix].astype(np.float64)
    # 30-min step: acc_mm = intensity_mmh * 0.5 h
    ann_sum = float(np.nansum(np.maximum(vals, 0.0)) * 0.5)
    rows.append({'year': yr, 'raster_mm': ann_sum})

df_r = pd.DataFrame(rows).set_index('year')

# Station data from calib CSV
df_s = pd.read_csv(CALIB_CSV)
df_s['datetime'] = pd.to_datetime(df_s['datetime'])
df_s['year'] = df_s['datetime'].dt.year
grp = df_s.groupby('year')[['P_station_mm', 'P_sat_mm', 'P_corrected_mm']].sum()

comp = df_r.join(grp, how='left')
comp['PBIAS_pct'] = (comp['raster_mm'] - comp['P_station_mm']) / comp['P_station_mm'] * 100

print(f"{'Year':>4}  {'Raster':>8}  {'Station':>8}  {'IMERG_raw':>10}  {'QM_csv':>8}  {'PBIAS%':>8}")
print("-" * 60)
for yr, row in comp[(comp.index >= 2001)].iterrows():
    print(f"{yr:>4}  {row.raster_mm:>8.1f}  {row.P_station_mm:>8.1f}  "
          f"{row.P_sat_mm:>10.1f}  {row.P_corrected_mm:>8.1f}  {row.PBIAS_pct:>+8.1f}%")

print("-" * 60)
valid = comp[(comp.index >= 2001) & comp['P_station_mm'].notna()]
print(f"Mean Raster: {valid.raster_mm.mean():.1f} mm/yr")
print(f"Mean Station: {valid.P_station_mm.mean():.1f} mm/yr")
print(f"Mean PBIAS: {valid.PBIAS_pct.mean():+.1f}%")
