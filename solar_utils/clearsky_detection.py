import pandas as pd
import numpy as np
import pvlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sklearn.model_selection import train_test_split

from sklearn.linear_model import LinearRegression, HuberRegressor

def ineichen(ghi_df, TL):
    """
    Calculate the clear sky irradiance using the Ineichen model.
    
    Parameters:
    ghi_df (pd.DataFrame): DataFrame containing GHI data with columns 'alt', 'cos_zna', and 'AM'.
    TL (float): Linke turbidity factor.
    
    Returns:
    np.ndarray: Calculated clear sky irradiance values.
    """
    I0 = 1366.1
    h = ghi_df['alt']
    f1 = np.exp(-h/8000)
    f2 = np.exp(-h/1250)
    a1 = (h*5.09e-5) + 0.868
    a2 = (h*3.92e-5) + 0.0387
    return a1*I0*ghi_df['cos_zna']*np.exp(-a2*ghi_df['AM']*(f1+f2*(TL-1)))

def cal_TL_nls_with_eval(clr_data_df, alt, test_size=0.2, random_state=42):
    """
    Calculate the Linke turbidity factor using non-linear least squares method with evaluation.

    Parameters:
    clr_data_df (pd.DataFrame): DataFrame containing clear sky data with columns 'Datetime', 'zna', 'I', and others.
    alt (float): Altitude of the location.
    test_size (float): Fraction of the data to be used for testing.
    random_state (int): Seed for the random number generator.

    Returns:
    tuple: TL_nls_array (np.ndarray), train_site_df (pd.DataFrame), test_site_df (pd.DataFrame)
    """
    h = alt
    TL_nls_list = []
    
    if eval:
        train_site_df = pd.DataFrame()
        test_site_df = pd.DataFrame()

    for month in range(1, 13):
        month_data = clr_data_df[clr_data_df['Datetime'].dt.month == month]
        month_data = month_data[month_data['zna']<85]
        month_data = month_data.dropna()
        month_data['alt'] = h
        
        if month_data.shape[0] >= 10:
            train_month_df, test_month_df = train_test_split(month_data, test_size=test_size, random_state=random_state)

            opt_TL, _ = curve_fit(ineichen, train_month_df, train_month_df['I'])
            TL_nls_list.append(opt_TL[0])

            train_month_df['month'] = month
            test_month_df['month'] = month

            train_site_df = pd.concat([train_site_df, train_month_df], axis=0, ignore_index=True)
            test_site_df = pd.concat([test_site_df, test_month_df], axis=0, ignore_index=True)
                
        elif 0 <= month_data.shape[0] < 10:
            train_month_df = month_data.copy()
            opt_TL, _ = curve_fit(ineichen, train_month_df, train_month_df['I'])
            TL_nls_list.append(opt_TL[0])

            train_site_df = pd.concat([train_site_df, train_month_df], axis=0, ignore_index=True)
        else:   
            TL_nls_list.append(np.nan)

    TL_nls_series = pd.Series(TL_nls_list)

    # Interpolate NaN values in the series using forward and backward fill
    TL_nls_series = TL_nls_series.interpolate(method='linear', limit_direction='both')

    # Convert the series to a numpy array
    TL_nls_array = np.array(TL_nls_series)

    # Reshape to ensure a 2D array with shape (n, 1)
    TL_nls_array = TL_nls_array.reshape(-1, 1)

    return TL_nls_array, train_site_df, test_site_df

def cal_TL_nls(clr_data_df, alt, graph_plot=False, path='TL_pics'):
    """
    Calculate the Linke turbidity factor using non-linear least squares method.

    Parameters:
    clr_data_df (pd.DataFrame): DataFrame containing clear sky data with columns 'Datetime', 'zna', 'I', and others.
    alt (float): Altitude of the location.
    graph_plot (bool): Whether to plot graphs of the estimated clear sky irradiance.
    path (str): Path to save the plots.

    Returns:
    np.ndarray: Calculated TL_nls_array.
    """
    h = alt
    TL_nls_list = []
    fig, axs = plt.subplots(3, 4, figsize=(45, 15))

    for month in range(1, 13):
        month_data = clr_data_df[clr_data_df['Datetime'].dt.month == month]
        month_data = month_data[month_data['zna']<85]
        month_data = month_data.dropna()
        month_data['alt'] = h
        
        if not month_data.empty:
            opt_TL, _ = curve_fit(ineichen, month_data, month_data['I'])
            TL_nls_list.append(opt_TL[0])

            if graph_plot:
                row, col = divmod(month-1, 4)
                iclr_hat = ineichen(month_data, TL=opt_TL)
    
                plot_df = month_data.copy()
                plot_df['est_Iclr'] = iclr_hat.values
                
                axs[row, col].plot(plot_df['I'].values, label='detected ghi', marker='o')
                axs[row, col].plot(plot_df['est_Iclr'].values, label=f'estimated Iclr with TL = {opt_TL[0]:.4f}', marker='x')
                axs[row, col].set_title(f"{month_data['site_name'].values[0]} Month {month}")
                axs[row, col].set_ylabel('Irradiance [W/m2]')
                axs[row, col].set_xlabel('Time step')
                axs[row, col].legend()
        else:
            TL_nls_list.append(np.nan)

    if graph_plot:
        fig.savefig(fname=f"{path}/{clr_data_df['site_name'].values[0]}_compare_estimateTL.png")
        
    TL_nls_series = pd.Series(TL_nls_list)

    # Interpolate NaN values in the series using forward and backward fill
    TL_nls_series = TL_nls_series.interpolate(method='linear', limit_direction='both')

    # Convert the series to a numpy array
    TL_nls_array = np.array(TL_nls_series)

    # Reshape to ensure a 2D array with shape (n, 1)
    TL_nls_array = TL_nls_array.reshape(-1, 1)

    return TL_nls_array

def detect_clr_portion_pvlib(site_df, window_length=45, max_diff=150, mean_diff=150, freq='15min', slope_dev=100, max_iterations=30):
    """
    Detect clear sky portions in GHI measurements using the pvlib library.

    Parameters:
    site_df (pd.DataFrame): DataFrame containing site data with columns 'Datetime', 'I', and 'Iclr'.
    window_length (int): Length of the window for clear sky detection.
    max_diff (float): Maximum difference threshold for clear sky detection.
    mean_diff (float): Mean difference threshold for clear sky detection.
    freq (str): Frequency of the time series data.
    slope_dev (float): Slope deviation threshold for clear sky detection.
    max_iterations (int): Maximum number of iterations for clear sky detection.

    Returns:
    tuple: detected_clr_ghi_df (pd.DataFrame), detected_partly_clr_day_df (pd.DataFrame), clr_date_list (list)
    
    Reference:
    Reno, M.J. and C.W. Hansen, “Identification of periods of clear sky irradiance in time series of GHI measurements” 
    Renewable Energy, v90, p. 520-531, 2016.
    """
    site_df['Datetime'] = pd.to_datetime(site_df['Datetime'])
    site_df = site_df.set_index("Datetime")

    all_date_list = sorted(list(set(site_df.index.date)))
    detected_clr_ghi_df = pd.DataFrame()
    detected_partly_clr_day_df = pd.DataFrame()
    clr_date_list = []

    for date in all_date_list:
        used_df = pd.DataFrame(index=pd.date_range(start=f"{date} 07:00", end=f"{date} 17:01", freq=freq))
        used_df.index.name = 'Datetime'
        used_df = used_df.merge(site_df, left_index=True, right_index=True, how='left')

        try:
            clr_samples = pvlib.clearsky.detect_clearsky(measured=used_df['I'],
                                                        clearsky=used_df['Iclr'],
                                                        max_diff=max_diff,
                                                        mean_diff=mean_diff,
                                                        window_length=window_length,
                                                        slope_dev=slope_dev,
                                                        max_iterations=max_iterations)
            
            if np.sum(clr_samples) != 0:
                used_df['is_clr_point'] = clr_samples
                clr_date_list.append(date)
                clr_portion = used_df[clr_samples]
                clr_portion = clr_portion.reset_index()

                detected_clr_ghi_df = pd.concat([detected_clr_ghi_df, clr_portion], axis=0, ignore_index=True)
                used_df = used_df.reset_index()
                detected_partly_clr_day_df = pd.concat([detected_partly_clr_day_df, used_df], axis=0, ignore_index=True)
        except Exception as e:
            print(f"{date} is terminated due to an error: {e}")
    
    return detected_clr_ghi_df, detected_partly_clr_day_df, clr_date_list
