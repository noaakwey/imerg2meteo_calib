import os
import glob
import pandas as pd
from tqdm import tqdm

from data_loader import load_meteo_station, load_imerg_data, pair_datasets
from qm_calibration import calibrate_station
from validation import aggregate_and_validate

def main():
    base_dir = r"d:\Cache\Yandex.Disk\РАЗРАБОТКА\code\imerg2meteo_calib"
    meteo_dir = os.path.join(base_dir, "data", "meteo", "срочные данные_осадки")
    imerg_dir = os.path.join(base_dir, "data", "imerg", "IMERG_STATIONS_3H")
    
    out_dir = os.path.join(base_dir, "output")
    out_dir_calib = os.path.join(out_dir, "calib")
    os.makedirs(out_dir_calib, exist_ok=True)
    
    print("Loading IMERG dataset (this may take a minute)...")
    imerg_all = load_imerg_data(imerg_dir)
    print(f"IMERG dataset loaded: {len(imerg_all)} non-NaN records.")
    
    meteo_files = glob.glob(os.path.join(meteo_dir, "*.csv"))
    print(f"Found {len(meteo_files)} meteo stations.")
    
    all_metrics = []
    
    # Process all stations
    limit_stations = len(meteo_files)
    
    for mf in tqdm(meteo_files[:limit_stations], desc="Processing stations"):
        st_name = os.path.splitext(os.path.basename(mf))[0]
        
        # 1. Load meteo
        try:
            meteo_df = load_meteo_station(mf)
        except Exception as e:
            print(f"Error loading {st_name}: {e}")
            continue
            
        if len(meteo_df) == 0:
            continue
            
        wmo_idx = meteo_df['wmo_index'].iloc[0]
        
        # 2. Filter IMERG to this station
        imerg_st = imerg_all[imerg_all['wmo_index'] == wmo_idx]
        if imerg_st.empty:
            continue
            
        # 3. Pair datasets
        paired_df = pair_datasets(meteo_df, imerg_st)
        
        # 4. Calibrate (wet-day QM)
        calibrated_df = calibrate_station(paired_df, train_start='2001-01-01', train_end='2015-12-31')
        
        # 5. Save calibration results
        out_csv = os.path.join(out_dir_calib, f"{st_name}_{wmo_idx}_calib.csv")
        calibrated_df.to_csv(out_csv, index=False)
        
        # 6. Evaluate metrics
        metrics = aggregate_and_validate(calibrated_df, date_start='2016-01-01', date_end='2021-12-31')
        
        # Flatten metrics for logging
        row = {
            'wmo_index': wmo_idx,
            'station_name': st_name,
            'daily_KGE_raw': metrics['daily']['KGE_raw'],
            'daily_KGE_corr': metrics['daily']['KGE_corr'],
            'daily_PBIAS_raw': metrics['daily']['PBIAS_raw'],
            'daily_PBIAS_corr': metrics['daily']['PBIAS_corr'],
            'monthly_KGE_raw': metrics['monthly']['KGE_raw'],
            'monthly_KGE_corr': metrics['monthly']['KGE_corr'],
            'monthly_PBIAS_raw': metrics['monthly']['PBIAS_raw'],
            'monthly_PBIAS_corr': metrics['monthly']['PBIAS_corr'],
        }
        all_metrics.append(row)
        
    # Save overall metrics
    if all_metrics:
        metrics_df = pd.DataFrame(all_metrics)
        metrics_csv = os.path.join(out_dir, "validation_metrics_test.csv")
        metrics_df.to_csv(metrics_csv, index=False)
        print(f"Metrics saved to {metrics_csv}")
        print("\nAverage metrics sample run:")
        print(metrics_df[['daily_KGE_raw', 'daily_KGE_corr', 'monthly_KGE_raw', 'monthly_KGE_corr']].mean())

if __name__ == '__main__':
    main()
