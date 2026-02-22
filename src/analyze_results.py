import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

def main():
    parser = argparse.ArgumentParser(description="Analyze calibration results.")
    parser.add_argument("--dataset", type=str, choices=['imerg', 'era5land'], default='imerg',
                        help="Choose which dataset to analyze (default: imerg)")
    args = parser.parse_args()
    
    base_dir = r"d:\Cache\Yandex.Disk\РАЗРАБОТКА\code\imerg2meteo_calib"
    metrics_file = os.path.join(base_dir, "output", f"calib_{args.dataset}", f"validation_metrics_{args.dataset}.csv")
        
    out_plots_dir = os.path.join(base_dir, "output", f"plots_{args.dataset}")
    os.makedirs(out_plots_dir, exist_ok=True)
    
    if not os.path.exists(metrics_file):
        print(f"Metrics file not found: {metrics_file}. Run main.py --dataset {args.dataset} first.")
        return
        
    df = pd.read_csv(metrics_file)
    
    # 1. Generate Summary Statistics Table
    cols = ['daily_KGE_raw', 'daily_KGE_corr', 'monthly_KGE_raw', 'monthly_KGE_corr', 
            'daily_PBIAS_raw', 'daily_PBIAS_corr', 'monthly_PBIAS_raw', 'monthly_PBIAS_corr']
            
    stats = df[cols].agg(['mean', 'median', 'std']).T
    stats.to_csv(os.path.join(out_plots_dir, "summary_statistics.csv"))
    print(f"Summary Statistics ({args.dataset.upper()}):")
    print(stats)
    
    # 2. Boxplots of KGE (Monthly and Daily)
    plt.figure(figsize=(10, 6))
    kge_data = pd.melt(df, value_vars=['daily_KGE_raw', 'daily_KGE_corr', 'monthly_KGE_raw', 'monthly_KGE_corr'], 
                       var_name='Metric', value_name='KGE')
    # Filter out extreme outliers for better visualization
    kge_data_filtered = kge_data[(kge_data['KGE'] > -2) & (kge_data['KGE'] < 1)]
    
    sns.boxplot(data=kge_data_filtered, x='Metric', y='KGE', palette='Set2')
    plt.title(f'Распределение метрики KGE до и после калибровки ({args.dataset.upper()} QM)')
    plt.ylabel('KGE (Kling-Gupta Efficiency)')
    plt.xticks(rotation=15)
    plt.axhline(0, color='red', linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(os.path.join(out_plots_dir, "kge_boxplot.png"), dpi=300)
    plt.close()
    
    # 3. Boxplots of PBIAS (Monthly and Daily)
    plt.figure(figsize=(10, 6))
    pbias_data = pd.melt(df, value_vars=['daily_PBIAS_raw', 'daily_PBIAS_corr', 'monthly_PBIAS_raw', 'monthly_PBIAS_corr'], 
                         var_name='Metric', value_name='PBIAS')
    pbias_data_filtered = pbias_data[(pbias_data['PBIAS'] > -100) & (pbias_data['PBIAS'] < 200)]
    
    sns.boxplot(data=pbias_data_filtered, x='Metric', y='PBIAS', palette='Set3')
    plt.title(f'Распределение ошибки PBIAS (%) до и после калибровки ({args.dataset.upper()} QM)')
    plt.ylabel('PBIAS (%)')
    plt.xticks(rotation=15)
    plt.axhline(0, color='red', linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(os.path.join(out_plots_dir, "pbias_boxplot.png"), dpi=300)
    plt.close()

    # 4. Scatter plot for the best station (Monthly)
    best_station = df.loc[df['monthly_KGE_corr'].idxmax()]
    best_wmo = best_station['wmo_index']
    best_name = best_station['station_name']
    
    calib_files = glob.glob(os.path.join(base_dir, "output", f"calib_{args.dataset}", f"*_{best_wmo}_calib.csv"))
    if calib_files:
        best_df = pd.read_csv(calib_files[0])
        best_df['datetime'] = pd.to_datetime(best_df['datetime'])
        best_df['year_month'] = best_df['datetime'].dt.to_period('M')
        
        # Monthly aggregates for the best station validation period (2016-2021)
        val_mask = (best_df['datetime'] >= '2016-01-01') & (best_df['datetime'] <= '2021-12-31')
        monthly_best = best_df[val_mask].groupby('year_month').agg({
            'P_station_mm': 'sum',
            'P_sat_mm': 'sum',
            'P_corrected_mm': 'sum'
        }).reset_index()
        
        plt.figure(figsize=(12, 5))
        
        plt.subplot(1, 2, 1)
        plt.scatter(monthly_best['P_station_mm'], monthly_best['P_sat_mm'], alpha=0.6, color='blue')
        plt.plot([0, monthly_best['P_station_mm'].max()], [0, monthly_best['P_station_mm'].max()], 'r--')
        plt.title(f'Сырые данные {args.dataset.upper()}\nСтанция {best_name}')
        plt.xlabel('Осадки по метеостанции, мм/месяц')
        plt.ylabel(f'Осадки {args.dataset.upper()} (сырые), мм/месяц')
        plt.grid(True, alpha=0.3)
        
        plt.subplot(1, 2, 2)
        plt.scatter(monthly_best['P_station_mm'], monthly_best['P_corrected_mm'], alpha=0.6, color='green')
        plt.plot([0, monthly_best['P_station_mm'].max()], [0, monthly_best['P_station_mm'].max()], 'r--')
        plt.title(f'Откалиброванные данные {args.dataset.upper()} (QM)\nСтанция {best_name}')
        plt.xlabel('Осадки по метеостанции, мм/месяц')
        plt.ylabel(f'Осадки {args.dataset.upper()} (калибр.), мм/месяц')
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(os.path.join(out_plots_dir, "scatter_comparison_best_station.png"), dpi=300)
        plt.close()

if __name__ == '__main__':
    main()
