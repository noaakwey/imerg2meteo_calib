"""
Legacy prototype script.

Production raster calibration pipeline is implemented in:
  src/apply_qm_to_rasters.py

Current IMERG main mode: v5_year_anchor.
"""

import os
import glob
import pandas as pd
import numpy as np
import xarray as xr
from scipy.spatial import cKDTree
from tqdm import tqdm

try:
    import rioxarray as rxr
except ImportError:
    print("rioxarray not installed. Please run: pip install rioxarray")

# Import our QM functions
from qm_calibration import get_season, fit_qm_transfer, apply_qm, fit_highres_qm

def load_station_metadata(meteo_dir):
    """
    Loads WMO index, Name, Longitude (X), and Latitude (Y) from Meteo CSVs.
    """
    csv_files = glob.glob(os.path.join(meteo_dir, "*.csv"))
    meta = []
    for f in csv_files:
        try:
            df = pd.read_csv(f, sep=';', encoding='cp866', nrows=1)
            meta.append({
                'wmo_index': df['Index'].iloc[0],
                'name': df['StationName'].iloc[0],
                'X': df['X'].iloc[0],
                'Y': df['Y'].iloc[0]
            })
        except Exception as e:
            pass
    return pd.DataFrame(meta)

def precalculate_qm_models(calib_dir, train_start='2001-01-01', train_end='2015-12-31', dataset='imerg'):
    """
    Reads the calibrated CSVs to extract training paired data and refits QM.
    Returns: dict[season][wmo_index] = (q_imerg, q_station, slope)
    """
    models = {'DJF': {}, 'MAM': {}, 'JJA': {}, 'SON': {}}
    calib_files = glob.glob(os.path.join(calib_dir, "*_calib.csv"))
    
    for f in tqdm(calib_files, desc="Precalculating QM models"):
        df = pd.read_csv(f)
        wmo = df['wmo_index'].iloc[0]
        
        df['datetime'] = pd.to_datetime(df['datetime'])
        train_mask = (df['datetime'] >= train_start) & (df['datetime'] <= train_end)
        train_df = df[train_mask]
        
        for season in models.keys():
            season_df = train_df[train_df['season'] == season]
            
            p_im_all = season_df['P_sat_mm'].values
            p_st_all = season_df['P_station_mm'].values
            
            if len(p_im_all) >= 20 and len(p_st_all) >= 20:
                if dataset == 'era5land':
                    q_sat, q_station, slope, p_th = fit_highres_qm(p_im_all, p_st_all, num_quantiles=1000)
                    if q_sat is not None:
                        models[season][wmo] = (q_sat, q_station, slope, p_th)
                else:
                    q_range = np.arange(0.01, 1.00, 0.01)
                    qm_params = fit_qm_transfer(p_im_all, p_st_all, q_range)
                    if qm_params[0] is not None:
                        models[season][wmo] = qm_params
                
    return models

def batch_process_rasters(tif_dir, out_dir, meteo_dir, calib_dir):
    os.makedirs(out_dir, exist_ok=True)
    
    print("1. Loading station metadata...")
    st_meta = load_station_metadata(meteo_dir)
    print(f"Loaded metadata for {len(st_meta)} stations.")
    
    print("2. Precalculating QM parameters...")
    qm_models = precalculate_qm_models(calib_dir)
    
    # List all TIFs
    tif_files = glob.glob(os.path.join(tif_dir, "*.tif"))
    if not tif_files:
        print("No TIF files found in directory.")
        return
        
    for tif_file in tif_files:
        fname = os.path.basename(tif_file)
        print(f"\nProcessing {fname}...")
        
        ds = rxr.open_rasterio(tif_file)
        # ds shape is usually (band, y, x)
        
        # 3. Build spatial Nearest Neighbor grid using KDTree
        # Get coordinates of the raster
        lons, lats = np.meshgrid(ds.x.values, ds.y.values)
        pixel_coords = np.column_stack([lons.ravel(), lats.ravel()])
        
        st_coords = st_meta[['X', 'Y']].values
        tree = cKDTree(st_coords)
        
        print("  - Building nearest station grid...")
        _, nearest_idx = tree.query(pixel_coords)
        nearest_wmo = st_meta['wmo_index'].values[nearest_idx]
        nearest_wmo_grid = nearest_wmo.reshape(lons.shape)
        
        # Parse bands into datetimes. Assuming bands are named P_YYYYMMDD_HHmm
        # We'll just generate a continuous 1-hour timeline assuming bands are consecutive.
        # Check if long_name contains the date.
        band_names = ds.attrs.get('long_name', ds.coords.get('band').values)
        if isinstance(band_names, str): 
            band_names = [band_names] # Fallback
            
        # Simplest robust way: if GEE exported it, bands might just be indices 1..N
        # Real scripts should parse timestamps from band names, but assuming they are sorted:
        # P_YYYYMMDD_HHmm structure -> try extracting from ds.band directly or provide start date.
        # For this template, we'll assume sequential 1-hour slots from the TIF filename context (e.g., 2008 Q3)
        # However, to be general, let's group dynamically by taking every 3 bands.
        
        n_bands = ds.shape[0]
        calibrated_bands = []
        
        print(f"  - Applying temporal downscaling ({n_bands} bands/hours)...")
        # Every 3 hours forms a specific 3-hour synoptic term window
        for i in tqdm(range(0, n_bands, 3), desc="3-hour windows", leave=False):
            # 1-hour arrays (up to 3)
            idx_end = min(i+3, n_bands)
            raw_1h = ds.values[i:idx_end, :, :] 
            
            # P_raw_3h
            raw_3h = np.sum(raw_1h, axis=0)
            
            # We need to determine the month to know the season.
            # In a full script, extract timestamp from band metadata. Here, we'll use a placeholder 'JJA' 
            # if we can't parse it. You must adapt date parsing based on your GEE exact output format.
            season = 'JJA' # Fallback. PLEASE REPLACE with actual parsing e.g. from band P_20080701_0000
            
            # Initialize calibrated 3h array
            calib_3h = np.zeros_like(raw_3h, dtype=float)
            
            # Apply QM per station
            unique_wmos = np.unique(nearest_wmo_grid)
            for wmo in unique_wmos:
                mask = (nearest_wmo_grid == wmo) & (raw_3h > 0)
                vals = raw_3h[mask]
                
                if len(vals) > 0:
                    params = qm_models.get(season, {}).get(wmo, None)
                    if params is not None:
                        # Both High-Res EQM and Standard QM emit a 4-tuple: (q_sat, q_station, slope, p_th)
                        q_im, q_st, slope, p_th = params
                        calib_3h[mask] = apply_qm(vals, q_im, q_st, slope, p_th)
                    else:
                        calib_3h[mask] = vals # fallback if no model
                        
            # Calculate weights and downscale
            weights = np.zeros_like(raw_1h, dtype=float)
            mask_gt0 = raw_3h > 0
            
            for j in range(raw_1h.shape[0]):
                weights[j, mask_gt0] = raw_1h[j, mask_gt0] / raw_3h[mask_gt0]
                
            calib_1h = calib_3h * weights
            calibrated_bands.append(calib_1h)
            
        print("  - Exporting GeoTIFF...")
        final_array = np.concatenate(calibrated_bands, axis=0)
        
        # Create output dataset preserving spatial coordinates
        out_ds = ds.copy(data=final_array)
        out_path = os.path.join(out_dir, fname.replace('.tif', '_qm_calibrated.tif'))
        
        # Ensure we write floats cleanly
        out_ds.rio.to_raster(out_path, compress='lzw')
        print(f"  - Saved to {out_path}")

if __name__ == '__main__':
    # Define paths
    BASE_DIR = r"d:\Cache\Yandex.Disk\РАЗРАБОТКА\code\imerg2meteo_calib"
    INPUT_RST_DIR = os.path.join(BASE_DIR, "data", "imerg", "geotiff_1h") # Folder with GEE TIFs
    OUT_RST_DIR = os.path.join(BASE_DIR, "output", "geotiff_calib")
    METEO_DIR = os.path.join(BASE_DIR, "data", "meteo", "срочные данные_осадки")
    CALIB_CSV_DIR = os.path.join(BASE_DIR, "output", "calib_imerg")
    
    # Uncomment to run when TIF folder is populated:
    # batch_process_rasters(INPUT_RST_DIR, OUT_RST_DIR, METEO_DIR, CALIB_CSV_DIR)
