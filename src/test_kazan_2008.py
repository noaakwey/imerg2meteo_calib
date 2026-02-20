import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
try:
    import rioxarray as rxr
except ImportError:
    print("rioxarray not installed. Please install it (pip install rioxarray).")
    import sys
    sys.exit(0)

# 1. Coordinates for Kazan (from meteo file)
# To be robust, let's read them dynamically, or we can just use 55.79, 49.11
lat_kazan = 55.80
lon_kazan = 49.18 # roughly Kazan station

def main():
    base_dir = r"d:\Cache\Yandex.Disk\РАЗРАБОТКА\code\imerg2meteo_calib"
    tif_file = r"C:\Users\artur\Downloads\IMERG_V07_P1h_mm_2008_Q3_permanent.tif"
    calib_file = os.path.join(base_dir, "output", "calib", "Казань_27595_calib.csv")
    out_dir = os.path.join(base_dir, "output", "kazan_test")
    os.makedirs(out_dir, exist_ok=True)
    
    if not os.path.exists(tif_file):
        print(f"TIF file not found: {tif_file}")
        return
        
    print("Loading 3-hour calibrated data for Kazan...")
    df_calib = pd.read_csv(calib_file)
    df_calib['datetime'] = pd.to_datetime(df_calib['datetime'])
    
    # Filter for July 2008 (Q3 is Jul, Aug, Sep)
    df_jul2008 = df_calib[(df_calib['datetime'] >= '2008-07-01') & (df_calib['datetime'] < '2008-08-01')].copy()
    
    print("Loading 1-hour TIF data...")
    ds = rxr.open_rasterio(tif_file)
    
    # Extract point data for Kazan
    lon = 49.18
    lat = 55.80
    point_data = ds.sel(x=lon, y=lat, method="nearest")
    
    # Extract times from band names (P_YYYYMMdd_HHmm)
    band_names = point_data.long_name if hasattr(point_data, 'long_name') else point_data.attrs.get('long_name', point_data.name)
    
    # Wait, rioxarray creates a 'band' dimension. We need to parse dates from the 'band' dimension coordinates or attributes.
    # Often, GEE exports have band names as a tuple or list in ds.attrs['long_name']
    
    # Let's see how xarray presents it
    band_coords = ds.coords['band'].values
    # GEE band names are usually stored somewhere else, but let's try reading the raster differently if needed
    
    # A cleaner way with rioxarray for GEE multi-band exports:
    # point_data.values is an array of shape (N_bands,)
    vals_1h = point_data.values
    
    # Assuming GEE bands are sequentially 1-hour intervals starting from 2008-07-01 00:00:00
    times_1h = pd.date_range(start='2008-07-01 00:00:00', periods=len(vals_1h), freq='1H')
    
    df_1h = pd.DataFrame({
        'datetime_1h': times_1h,
        'P_raw_1h': vals_1h
    })
    
    # Group 1-hour data into 3-hour chunks to align with 3-hour calibration (00, 03, 06...)
    # IMERG 3-hour ending at 03:00 covers 00:00, 01:00, 02:00? Actually IMERG 3h slots are 
    # named by the *start* or *end* of the period. The user mentioned "P3H".
    # Usually 00, 03, 06.
    
    # Let's map each 1-hour slot to its corresponding 3-hour window
    # e.g., 00:00, 01:00, 02:00 belongs to the 03:00 timeframe.
    df_1h['window_3h'] = df_1h['datetime_1h'].dt.ceil('3H')
    
    # Merge with 3-hour calibrated data
    df_merged = pd.merge(df_1h, df_jul2008[['datetime', 'P_station_mm', 'P_imerg_mm', 'P_corrected_mm']], 
                         left_on='window_3h', right_on='datetime', how='left')
                         
    # Calculate 1-hour calibrated values
    # P_calib_1h = P_corr_3h * (P_raw_1h / P_raw_3h)
    
    # First, calculate P_raw_3h from the 1-hour data as the sum within the window
    # Wait, the df_merged already has P_imerg_mm (which is the 3h sum from the calibrated file)
    # But let's calculate the sum of 1-hour raw to be safe as weights
    window_sums = df_1h.groupby('window_3h')['P_raw_1h'].sum().rename('P_raw_1h_sum')
    df_merged = df_merged.merge(window_sums, on='window_3h', how='left')
    
    # Avoid div zero
    df_merged['weight'] = 0.0
    mask_gt0 = df_merged['P_raw_1h_sum'] > 0
    df_merged.loc[mask_gt0, 'weight'] = df_merged.loc[mask_gt0, 'P_raw_1h'] / df_merged.loc[mask_gt0, 'P_raw_1h_sum']
    
    # Where 1-hour sum is zero but 3-hour corrected is not (rare), distribute evenly
    mask_eq0_corr_gt0 = (df_merged['P_raw_1h_sum'] == 0) & (df_merged['P_corrected_mm'] > 0)
    df_merged.loc[mask_eq0_corr_gt0, 'weight'] = 1.0 / 3.0
    
    df_merged['P_calib_1h'] = df_merged['P_corrected_mm'] * df_merged['weight']
    
    # Analyze July 21st, 2008 extreme event
    jul21 = df_merged[(df_merged['datetime_1h'] >= '2008-07-21 00:00:00') & 
                      (df_merged['datetime_1h'] < '2008-07-22 00:00:00')]
                      
    print("\n--- Июль 21, 2008 (г. Казань) ---")
    print(f"Сырой IMERG сумма (1h aggregated) за день: {jul21['P_raw_1h'].sum():.2f} мм")
    print(f"Сырой IMERG сумма (3h) за день: {jul21.drop_duplicates('window_3h')['P_imerg_mm'].sum():.2f} мм")
    print(f"Откалиброванный IMERG сумма (QM) за день: {jul21['P_calib_1h'].sum():.2f} мм")
    print(f"Фактические осадки по станции: {jul21.drop_duplicates('window_3h')['P_station_mm'].sum():.2f} мм")
    
    print("\nМаксимальная часовая интенсивность (I30 proxy / I60) 21 июля:")
    max_raw = jul21['P_raw_1h'].max()
    max_calib = jul21['P_calib_1h'].max()
    print(f"Сырая: {max_raw:.2f} мм/час")
    print(f"Откалиброванная: {max_calib:.2f} мм/час")
    
    # Monthly Total
    print(f"\n--- Июль 2008 в целом (г. Казань) ---")
    print(f"Сырой IMERG сумма за июль: {df_merged['P_raw_1h'].sum():.2f} мм")
    print(f"Откалиброванный IMERG (QM): {df_merged['P_calib_1h'].sum():.2f} мм")
    station_monthly = df_jul2008['P_station_mm'].sum()
    print(f"Фактические по станции: {station_monthly:.2f} мм")

if __name__ == '__main__':
    main()
