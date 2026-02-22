import numpy as np
import pandas as pd
import scipy.stats as stats

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

def fit_highres_qm(p_sat_all, p_st_all, num_quantiles=1000):
    """
    Fits High-Resolution Empirical Quantile Mapping (EQM) for Big Data (e.g. ERA5).
    Uses a dense array of quantiles to map the shape without parametric smoothing.
    """
    # 1. Frequency Adaptation (Calculate PoP)
    pop_st = np.sum(p_st_all > 0) / len(p_st_all) if len(p_st_all) > 0 else 0
    
    if pop_st == 0:
        return None, None, 1.0, 0.0
        
    p_th = np.quantile(p_sat_all, 1.0 - pop_st)
    if p_th < 0:
        p_th = 0
        
    p_sat_wet = p_sat_all[p_sat_all > p_th]
    p_st_wet = p_st_all[p_st_all > 0]
    
    if len(p_sat_wet) < 20 or len(p_st_wet) < 20:
        return None, None, 1.0, p_th
        
    # Create high-res linear space (e.g., 0.001 to 0.999)
    # Using linspace avoids hard 0 and 1 which can cause infinity in parametric, though eqm handles it fine
    quantiles_list = np.linspace(0.001, 0.999, num_quantiles)
    
    q_sat = compute_quantiles(p_sat_wet, quantiles_list)
    q_station = compute_quantiles(p_st_wet, quantiles_list)
    
    # Linear extrapolation for the extreme right tail (beyond 99.9%)
    if (q_sat[-1] - q_sat[-5]) > 0:
        dy = q_station[-1] - q_station[-5]
        dx = q_sat[-1] - q_sat[-5]
        slope_extrap = dy / dx if dx > 0 else 1.0
    else:
        slope_extrap = 1.0
        
    return q_sat, q_station, slope_extrap, p_th

def fit_qm_transfer(p_sat_all, p_st_all, quantiles_list=np.arange(0.01, 1.00, 0.01)):
    """
    Fits standard Wet-Day QM transfer function between Satellite and station (For IMERG).
    """
    p_th = 0.0
    p_sat_wet = p_sat_all[p_sat_all > 0]
    p_st_wet = p_st_all[p_st_all > 0]
    
    if len(p_sat_wet) < 5 or len(p_st_wet) < 5:
        return None, None, 1.0, p_th
        
    # Calculate quantiles
    q_sat = compute_quantiles(p_sat_wet, quantiles_list)
    q_station = compute_quantiles(p_st_wet, quantiles_list)
    
    # Linear extrapolation for values above 99th percentile
    if len(p_sat_wet) > 20 and len(p_st_wet) > 20 and (q_sat[-1] - q_sat[-5]) > 0:
        dy = q_station[-1] - q_station[-5]
        dx = q_sat[-1] - q_sat[-5]
        slope_extrap = dy / dx if dx > 0 else 1.0
    else:
        slope_extrap = 1.0
        
    return q_sat, q_station, slope_extrap, p_th

def apply_qm(val, q_sat, q_station, slope_extrap, p_th):
    """
    Applies the fitted QM transformation with precipitation thresholding (p_th).
    """
    if isinstance(val, (int, float)):
        return _apply_qm_scalar(val, q_sat, q_station, slope_extrap, p_th)
    elif isinstance(val, pd.Series):
        val_arr = val.values
        res = np.zeros_like(val_arr, dtype=float)
        
        # Values strictly greater than threshold get mapped
        pos_mask = val_arr > p_th
        res[pos_mask] = [_apply_qm_scalar(v, q_sat, q_station, slope_extrap, p_th) for v in val_arr[pos_mask]]
        return res
    else:
        # np array
        res = np.zeros_like(val, dtype=float)
        pos_mask = val > p_th
        res[pos_mask] = [_apply_qm_scalar(v, q_sat, q_station, slope_extrap, p_th) for v in val[pos_mask]]
        return res

def _apply_qm_scalar(v, q_sat, q_station, slope_extrap, p_th):
    if v <= p_th:
        return 0.0
    if v <= q_sat[0]:
        # Interp between threshold and 1st quantile
        # Map values between p_th and q_sat[0] linearly to 0 -> q_station[0]
        dx = q_sat[0] - p_th
        if dx > 0:
            return q_station[0] * ((v - p_th) / dx)
        else:
            return q_station[0]
            
    if v >= q_sat[-1]:
        # Extrapolate beyond 99th
        return q_station[-1] + slope_extrap * (v - q_sat[-1])
    
    # Interpolate using binary search
    idx = np.searchsorted(q_sat, v)
    if idx == 0: return q_station[0] # Fallback
    
    x0, x1 = q_sat[idx-1], q_sat[idx]
    y0, y1 = q_station[idx-1], q_station[idx]
    
    if x1 == x0:
        return y0
    
    return y0 + (y1 - y0) * ((v - x0) / (x1 - x0))

def calibrate_station(paired_df, train_start='2001-01-01', train_end='2015-12-31', dataset='imerg'):
    """
    Full pipeline for a single station. Uses Parametric Gamma (ERA5) or Empirical QM (IMERG).
    Adds a `P_corrected_mm` column.
    """
    df = paired_df.copy()
    if 'season' not in df.columns:
        df['season'] = df['datetime'].dt.month.apply(get_season)
        
    train_mask = (df['datetime'] >= train_start) & (df['datetime'] <= train_end)
    
    df['P_corrected_mm'] = df['P_sat_mm'] # default fallback
    
    for season in ['DJF', 'MAM', 'JJA', 'SON']:
        train_season_df = df[train_mask & (df['season'] == season)]
        p_sat_all = train_season_df['P_sat_mm'].values
        p_st_all = train_season_df['P_station_mm'].values
        
        if len(p_sat_all) < 20 or len(p_st_all) < 20:
            continue
            
        season_mask = df['season'] == season
        
        if dataset == 'era5land':
            # ============================================================
            # Full-Distribution QM (маппинг всего CDF целиком, не только wet-day)
            # ============================================================
            # Используем ВСЕ значения (включая нули и морось), чтобы не терять
            # объём через агрессивное обнуление порога p_th.
            
            num_q = 1000
            quantiles_list = np.linspace(0.0, 1.0, num_q + 2)[1:-1]  # исключаем 0.0 и 1.0
            
            q_sat = np.quantile(p_sat_all, quantiles_list)
            q_st = np.quantile(p_st_all, quantiles_list)
            
            # Наклон для экстраполяции хвоста (выше 99.9-го квантиля)
            if (q_sat[-1] - q_sat[-5]) > 0:
                slope_extrap = (q_st[-1] - q_st[-5]) / (q_sat[-1] - q_sat[-5])
            else:
                slope_extrap = 1.0
            
            # Применяем маппинг ко всему сезону
            corrected = apply_qm(df.loc[season_mask, 'P_sat_mm'], q_sat, q_st, slope_extrap, p_th=0.0)
            corrected = np.maximum(corrected, 0.0)  # физический порог: осадки >= 0
            
            # Volume Scaling по тренировочному периоду
            train_corr = apply_qm(train_season_df['P_sat_mm'], q_sat, q_st, slope_extrap, p_th=0.0)
            train_corr = np.maximum(train_corr, 0.0)
            mean_st = np.mean(p_st_all)
            mean_corr = np.mean(train_corr)
            
            if mean_corr > 0:
                volume_factor = mean_st / mean_corr
                volume_factor = np.clip(volume_factor, 0.5, 2.0)
            else:
                volume_factor = 1.0
            
            corrected = corrected * volume_factor
            df.loc[season_mask, 'P_corrected_mm'] = corrected
            
        else:
            # Full-Distribution QM для IMERG (1000 квантилей, весь CDF + Volume Scaling)
            num_q = 1000
            quantiles_list = np.linspace(0.0, 1.0, num_q + 2)[1:-1]

            q_sat = np.quantile(p_sat_all, quantiles_list)
            q_st  = np.quantile(p_st_all,  quantiles_list)

            if (q_sat[-1] - q_sat[-5]) > 0:
                slope_extrap = (q_st[-1] - q_st[-5]) / (q_sat[-1] - q_sat[-5])
            else:
                slope_extrap = 1.0

            corrected = apply_qm(df.loc[season_mask, 'P_sat_mm'],
                                 q_sat, q_st, slope_extrap, p_th=0.0)
            corrected = np.maximum(corrected, 0.0)

            # Volume Scaling по тренировочному периоду
            train_corr = apply_qm(train_season_df['P_sat_mm'],
                                  q_sat, q_st, slope_extrap, p_th=0.0)
            train_corr = np.maximum(train_corr, 0.0)
            mean_st   = np.mean(p_st_all)
            mean_corr = np.mean(train_corr)

            if mean_corr > 0:
                vf = np.clip(mean_st / mean_corr, 0.5, 2.0)
            else:
                vf = 1.0

            df.loc[season_mask, 'P_corrected_mm'] = corrected * vf
        
    return df


def calibrate_station_moving_window(paired_df, half_window=15, val_start=None, val_end=None):
    """
    Moving Window Full-Distribution QM для ERA5-Land.
    Для каждого года Y строит QM на данных из окна [Y-half_window, Y+half_window].
    При кросс-валидации (val_start/val_end заданы) данные валидационного периода
    исключаются из обучающего окна.
    """
    df = paired_df.copy()
    if 'season' not in df.columns:
        df['season'] = df['datetime'].dt.month.apply(get_season)
    
    df['P_corrected_mm'] = df['P_sat_mm']  # fallback
    df['year'] = df['datetime'].dt.year
    
    all_years = sorted(df['year'].unique())
    min_year, max_year = all_years[0], all_years[-1]
    
    num_q = 1000
    quantiles_list = np.linspace(0.0, 1.0, num_q + 2)[1:-1]
    
    for target_year in all_years:
        # Определяем границы окна
        win_start = max(min_year, target_year - half_window)
        win_end = min(max_year, target_year + half_window)
        
        # Маска обучающего окна
        window_mask = (df['year'] >= win_start) & (df['year'] <= win_end)
        
        # При кросс-валидации: исключаем валидационный период из обучения
        if val_start is not None and val_end is not None:
            val_s = int(val_start[:4])
            val_e = int(val_end[:4])
            window_mask = window_mask & ~((df['year'] >= val_s) & (df['year'] <= val_e))
        
        target_mask = df['year'] == target_year
        
        for season in ['DJF', 'MAM', 'JJA', 'SON']:
            train_data = df[window_mask & (df['season'] == season)]
            p_sat_win = train_data['P_sat_mm'].values
            p_st_win = train_data['P_station_mm'].values
            
            if len(p_sat_win) < 30 or len(p_st_win) < 30:
                continue
            
            target_season_mask = target_mask & (df['season'] == season)
            if target_season_mask.sum() == 0:
                continue
            
            # Full-Distribution QM (1000 квантилей)
            q_sat = np.quantile(p_sat_win, quantiles_list)
            q_st = np.quantile(p_st_win, quantiles_list)
            
            # Экстраполяция хвоста
            if (q_sat[-1] - q_sat[-5]) > 0:
                slope_extrap = (q_st[-1] - q_st[-5]) / (q_sat[-1] - q_sat[-5])
            else:
                slope_extrap = 1.0
            
            # Применяем к целевому году
            corrected = apply_qm(df.loc[target_season_mask, 'P_sat_mm'],
                                 q_sat, q_st, slope_extrap, p_th=0.0)
            corrected = np.maximum(corrected, 0.0)
            
            # Volume Scaling по окну обучения
            train_corr = apply_qm(train_data['P_sat_mm'],
                                  q_sat, q_st, slope_extrap, p_th=0.0)
            train_corr = np.maximum(train_corr, 0.0)
            mean_st = np.mean(p_st_win)
            mean_corr = np.mean(train_corr)
            
            if mean_corr > 0:
                vf = np.clip(mean_st / mean_corr, 0.5, 2.0)
            else:
                vf = 1.0
            
            df.loc[target_season_mask, 'P_corrected_mm'] = corrected * vf
    
    df.drop(columns=['year'], inplace=True)
    return df
