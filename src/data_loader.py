import pandas as pd
import numpy as np
import os
import glob
from pathlib import Path

def load_meteo_station(filepath):
    """
    Loads meteo station data from CSV.
    Retains only 3-hour accumulation intervals for direct comparison with IMERG 
    (hours: 3, 6, 15, 18).
    Keeps 00 and 12 UTC for daily/monthly controls.
    """
    df = pd.read_csv(filepath, sep=';', encoding='cp866')
    
    # Required columns: Index, Year, Month, Day, Hour, Precipitation, StationName, X, Y, H
    df['datetime'] = pd.to_datetime(
        dict(year=df['Year'], month=df['Month'], day=df['Day'], hour=df['Hour'])
    )
    
    df = df.rename(columns={
        'Index': 'wmo_index',
        'Precipitation': 'P_station_mm',
        'StationName': 'station_name'
    })
    
    # Extract only necessary columns
    df = df[['wmo_index', 'station_name', 'datetime', 'Hour', 'P_station_mm', 'X', 'Y', 'H']]
    return df

def load_satellite_data(folder_path, file_pattern):
    """
    Loads all satellite/reanalysis CSVs from the provided folder.
    Concatenates them into a single dataframe.
    """
    all_files = glob.glob(os.path.join(folder_path, file_pattern))
    
    df_list = []
    for f in all_files:
        _df = pd.read_csv(f)
        df_list.append(_df)
        
    sat_df = pd.concat(df_list, ignore_index=True)
    
    sat_df['datetime'] = pd.to_datetime(sat_df['time_utc'])
    sat_df = sat_df.rename(columns={
        'Index': 'wmo_index',
        'Name': 'station_name',
        'P_mm_3h': 'P_sat_mm'
    })
    
    # Drop NaNs: missing slots should be excluded from pairwise comparison
    sat_df = sat_df.dropna(subset=['P_sat_mm'])
    
    return sat_df[['wmo_index', 'datetime', 'P_sat_mm']]

def pair_datasets(meteo_df, sat_df):
    """
    Pairs satellite data to station data on wmo_index and datetime.
    Filters to only 3-hour terms (03, 06, 15, 18 UTC).
    For these terms, if Satellite has a record but Meteo is missing, Meteo gets P = 0.
    """
    # Find which 3-hour hours are supported by this specific station
    # (some only support 3 and 15, others also 6 and 18)
    supported_3h_hours = set(meteo_df['Hour'].unique()).intersection({3, 6, 15, 18})
    
    # Filter Satellite to only those supported hours
    sat_filtered = sat_df[sat_df['datetime'].dt.hour.isin(supported_3h_hours)].copy()
    
    # Filter meteo to these supported hours
    meteo_filtered = meteo_df[meteo_df['Hour'].isin(supported_3h_hours)][['wmo_index', 'datetime', 'P_station_mm']]
    
    # Merge left on Satellite so missing meteo records become NaNs
    merged = pd.merge(
        sat_filtered, 
        meteo_filtered, 
        on=['wmo_index', 'datetime'], 
        how='left'
    )
    
    # Missing meteo records indicate P_station = 0
    merged['P_station_mm'] = merged['P_station_mm'].fillna(0.0)
    
    return merged
