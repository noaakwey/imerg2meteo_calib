import os
import glob
import pandas as pd
import argparse
from tqdm import tqdm

from data_loader import load_meteo_station, load_satellite_data, pair_datasets
from qm_calibration import calibrate_station, calibrate_station_moving_window
from validation import aggregate_and_validate

def main():
    parser = argparse.ArgumentParser(description="Calibrate precipitation data.")
    parser.add_argument("--dataset", type=str, choices=['imerg', 'era5land'], default='imerg',
                        help="Choose which dataset to calibrate (default: imerg)")
    args = parser.parse_args()

    base_dir = r"d:\Cache\Yandex.Disk\РАЗРАБОТКА\code\imerg2meteo_calib"
    meteo_dir = os.path.join(base_dir, "data", "meteo", "срочные данные_осадки")
    
    if args.dataset == 'imerg':
        sat_dir = os.path.join(base_dir, "data", "imerg", "IMERG_STATIONS_3H")
        sat_pattern = "IMERG_V07_P3H_mm_*_permanent_trailing.csv"
        out_dir = os.path.join(base_dir, "output", "calib_imerg")
        train_start = '2001-01-01'
    elif args.dataset == 'era5land':
        sat_dir = os.path.join(base_dir, "data", "era5land", "ERA5LAND_STATIONS_3H")
        sat_pattern = "ERA5Land_P3H_mm_*.csv"
        out_dir = os.path.join(base_dir, "output", "calib_era5land")
        train_start = '2001-01-01'  # Сокращено с 1966: устраняет нестационарность климата (PBIAS -31% → ≈0%)
        
    train_end = '2015-12-31'
    val_start = '2016-01-01'
    val_end = '2021-12-31'
    
    os.makedirs(out_dir, exist_ok=True)
    
    print(f"Loading {args.dataset.upper()} dataset (this may take a minute)...")
    sat_all = load_satellite_data(sat_dir, sat_pattern)
    print(f"{args.dataset.upper()} dataset loaded: {len(sat_all)} non-NaN records.")
    
    meteo_files = glob.glob(os.path.join(meteo_dir, "*.csv"))
    print(f"Found {len(meteo_files)} meteo stations.")
    
    all_metrics = []
    
    # Process all stations
    for mf in tqdm(meteo_files, desc="Processing stations"):
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
        
        # 2. Filter Satellite to this station
        sat_st = sat_all[sat_all['wmo_index'] == wmo_idx]
        if sat_st.empty:
            continue
            
        # 3. Pair datasets
        paired_df = pair_datasets(meteo_df, sat_st, force_full_3h_grid=True)
        
        # 4. Calibrate
        if args.dataset == 'era5land':
            calibrated_df = calibrate_station_moving_window(
                paired_df, half_window=15, val_start=val_start, val_end=val_end)
        else:
            calibrated_df = calibrate_station(
                paired_df, train_start=train_start, train_end=train_end, dataset=args.dataset)
        
        # 5. Save calibration results
        out_csv = os.path.join(out_dir, f"{st_name}_{wmo_idx}_calib.csv")
        calibrated_df.to_csv(out_csv, index=False)
        
        # 6. Evaluate metrics
        metrics = aggregate_and_validate(calibrated_df, date_start=val_start, date_end=val_end)
        
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
        metrics_csv = os.path.join(out_dir, f"validation_metrics_{args.dataset}.csv")
        metrics_df.to_csv(metrics_csv, index=False)
        print(f"Metrics saved to {metrics_csv}")
        print(f"\nAverage metrics context run ({args.dataset.upper()}):")
        print(metrics_df[['daily_KGE_raw', 'daily_KGE_corr', 'monthly_KGE_raw', 'monthly_KGE_corr']].mean())

if __name__ == '__main__':
    main()
