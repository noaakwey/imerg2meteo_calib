import os
import pandas as pd
from data_loader import load_meteo_station, load_imerg_data, pair_datasets

imerg_path = r"d:\Cache\Yandex.Disk\РАЗРАБОТКА\code\imerg2meteo_calib\data\imerg\IMERG_STATIONS_3H"
meteo_file = r"d:\Cache\Yandex.Disk\РАЗРАБОТКА\code\imerg2meteo_calib\data\meteo\срочные данные_осадки\Акбулак.csv"

# Load one station
df_meteo = load_meteo_station(meteo_file)
print(f"Meteo: {len(df_meteo)} nonzero records loaded. Supported hours: {df_meteo['Hour'].unique()}")

# Loading only one year of IMERG for quick test to avoid long loading
test_imerg_file = os.path.join(imerg_path, "IMERG_V07_P3H_mm_2001_permanent_trailing.csv")
imerg_df = pd.read_csv(test_imerg_file)
imerg_df['datetime'] = pd.to_datetime(imerg_df['time_utc'])
imerg_df = imerg_df.rename(columns={'Index': 'wmo_index', 'Name': 'station_name', 'P_mm_3h': 'P_imerg_mm'})
imerg_df = imerg_df.dropna(subset=['P_imerg_mm'])
imerg_df = imerg_df[['wmo_index', 'datetime', 'P_imerg_mm']]
print(f"IMERG 2001: {len(imerg_df)} records loaded.")

# Filter IMERG to this station
imerg_df_station = imerg_df[imerg_df['wmo_index'] == df_meteo['wmo_index'].iloc[0]]

paired = pair_datasets(df_meteo, imerg_df_station)
print(f"Paired dataset: {len(paired)} rows.")
print("Sample pairs (P_station > 0):")
print(paired[paired['P_station_mm'] > 0].head())
print("\nSample pairs (P_station == 0):")
print(paired[paired['P_station_mm'] == 0].head())

# Check metrics
non_zero_meteo_count = (paired['P_station_mm'] > 0).sum()
non_zero_imerg_count = (paired['P_imerg_mm'] > 0).sum()
print(f"Non-zero meteo in 2001 for chosen hours: {non_zero_meteo_count}")
print(f"Non-zero IMERG in 2001 for chosen hours: {non_zero_imerg_count}")
