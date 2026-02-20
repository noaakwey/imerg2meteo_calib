import numpy as np
import pandas as pd

def get_season(month):
    if month in [12, 1, 2]:
        return 'DJF'
    elif month in [3, 4, 5]:
        return 'MAM'
    elif month in [6, 7, 8]:
        return 'JJA'
    else:
        return 'SON'

def compute_quantiles(series, q_range):
    """
    Computes exact quantiles for a series.
    series: pandas Series of non-zero values.
    q_range: array of probabilities [0.01, 0.02, ..., 0.99]
    """
    if len(series) == 0:
        return np.zeros_like(q_range)
    return np.quantile(series, q_range)

def fit_qm_transfer(p_imerg, p_station, quantiles_list=np.arange(0.01, 1.00, 0.01)):
    """
    Fits wet-day QM transfer function between IMERG and station for a specific season.
    p_imerg: array of IMERG non-zero values in training period
    p_station: array of station non-zero values in training period
    """
    # Calculate quantiles
    q_imerg = compute_quantiles(p_imerg, quantiles_list)
    q_station = compute_quantiles(p_station, quantiles_list)
    
    # Linear extrapolation for values above 99th percentile
    # Calculate slope based on the last few percentiles (e.g. 95th to 99th)
    if len(p_imerg) > 20 and len(p_station) > 20 and (q_imerg[-1] - q_imerg[-5]) > 0:
        dy = q_station[-1] - q_station[-5]
        dx = q_imerg[-1] - q_imerg[-5]
        slope_extrap = dy / dx if dx > 0 else 1.0
    else:
        slope_extrap = 1.0
        
    return q_imerg, q_station, slope_extrap

def apply_qm(val, q_imerg, q_station, slope_extrap):
    """
    Applies the fitted QM transformation to a single value or array of values.
    """
    if isinstance(val, (int, float)):
        return _apply_qm_scalar(val, q_imerg, q_station, slope_extrap)
    elif isinstance(val, pd.Series):
        val_arr = val.values
        res = np.zeros_like(val_arr, dtype=float)
        
        # Zero preservation
        pos_mask = val_arr > 0
        res[pos_mask] = [_apply_qm_scalar(v, q_imerg, q_station, slope_extrap) for v in val_arr[pos_mask]]
        return res
    else:
        # np array
        res = np.zeros_like(val, dtype=float)
        pos_mask = val > 0
        res[pos_mask] = [_apply_qm_scalar(v, q_imerg, q_station, slope_extrap) for v in val[pos_mask]]
        return res

def _apply_qm_scalar(v, q_imerg, q_station, slope_extrap):
    if v <= 0:
        return 0.0
    if v <= q_imerg[0]:
        # Interp between 0 and 1st quantile
        if q_imerg[0] > 0:
            return v * (q_station[0] / q_imerg[0])
        else:
            return v
    if v >= q_imerg[-1]:
        # Extrapolate beyond 99th
        return q_station[-1] + slope_extrap * (v - q_imerg[-1])
    
    # Interpolate
    idx = np.searchsorted(q_imerg, v)
    if idx == 0: return q_station[0] # Fallback
    
    x0, x1 = q_imerg[idx-1], q_imerg[idx]
    y0, y1 = q_station[idx-1], q_station[idx]
    
    if x1 == x0:
        return y0
    
    return y0 + (y1 - y0) * ((v - x0) / (x1 - x0))

def calibrate_station(paired_df, train_start='2001-01-01', train_end='2015-12-31'):
    """
    Full QM pipeline for a single station dataset containing `P_imerg_mm` and `P_station_mm`.
    Adds a `P_corrected_mm` column.
    """
    df = paired_df.copy()
    if 'season' not in df.columns:
        df['season'] = df['datetime'].dt.month.apply(get_season)
        
    train_mask = (df['datetime'] >= train_start) & (df['datetime'] <= train_end)
    
    df['P_corrected_mm'] = df['P_imerg_mm'] # default fallback
    
    for season in ['DJF', 'MAM', 'JJA', 'SON']:
        train_season_df = df[train_mask & (df['season'] == season)]
        
        # Wet-day QM: > 0 for both IMERG and Station separately
        p_im = train_season_df[train_season_df['P_imerg_mm'] > 0]['P_imerg_mm'].values
        p_st = train_season_df[train_season_df['P_station_mm'] > 0]['P_station_mm'].values
        
        if len(p_im) < 5 or len(p_st) < 5:
            # Not enough data to fit quantiles properly
            continue
            
        quantiles = np.arange(0.01, 1.00, 0.01)
        q_imerg, q_station, slope = fit_qm_transfer(p_im, p_st, quantiles)
        
        # Apply to all periods (train and valid)
        season_mask = df['season'] == season
        
        corrected = apply_qm(df.loc[season_mask, 'P_imerg_mm'], q_imerg, q_station, slope)
        df.loc[season_mask, 'P_corrected_mm'] = corrected
        
    return df
