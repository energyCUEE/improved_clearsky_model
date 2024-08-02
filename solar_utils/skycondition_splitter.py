import pandas as pd
import numpy as np

def split_sky_condition(site_df, method='ROC', min_threshold_k=0.6, clr_threshold_condition=0.8, partly_threshold_k=0.6):
    '''
    Splits the sky condition into clear, partly cloudy, and cloudy categories based on the provided method and thresholds.

    Parameters:
    site_df (pd.DataFrame): DataFrame containing site data with columns 'Datetime', 'cos_zna', 'I', 'Iclr', and 'k_bar'.
    method (str): Method to determine clear sky condition ('ROC' or 'k_bar'). Default is 'ROC'.
    min_threshold_k (float): Minimum threshold for the 'k' value to consider a day as clear of 'ROC' method. Default is 0.6.
    clr_threshold_condition (float): Threshold for 'k_bar' to consider a day as clear when method is 'k_bar'. Default is 0.8.
    partly_threshold_k (float): Threshold for 'k_bar' to consider a day as partly cloudy. Default is 0.6.

    Returns:
    tuple: Contains four lists:
        - clr_date_list (list): List of dates considered clear.
        - partly_cloudy_date_list (list): List of dates considered partly cloudy.
        - cloudy_date_list (list): List of dates considered cloudy.
        - average_k_list (list): List of average 'k' values for each date.
    '''

    # Parse to Datetime type then set as index
    site_df['Datetime'] = pd.to_datetime(site_df['Datetime'])
    site_df.set_index('Datetime', inplace=True)
    
    # Filter only day time
    site_df = site_df.between_time("07:00", "17:00") 
    
    # Consider only the time that degree of ZNA < 85
    site_df.loc[:,'zna'] = np.arccos(site_df['cos_zna']) * 180 / np.pi
    site_df = site_df[site_df['zna'] <= 85]

    site_df['k'] = site_df['I'] / site_df['Iclr']

    # Create list variables to stroing date informations 
    all_date_list = sorted(list(set(site_df.index.date)))
    clr_date_list = []
    partly_cloudy_date_list = []
    cloudy_date_list = []
    average_k_list = []

    for date in all_date_list:
        date_df = site_df[site_df.index.date == date].copy()  # Create a copy to avoid SettingWithCopyWarnings
        k_bar = date_df['k_bar'].iloc[0]
        date_df['is_increasing'] = (date_df['I'].diff() > 0).astype(int)

        date_df['is_concave_point'] = date_df['is_increasing'].diff()
        date_df['is_concave_point'] = date_df['is_concave_point'] == -1

        if method == 'ROC':
            clr_condition = (np.sum(date_df['is_concave_point']) == 1) and (np.sum(date_df['k'] <= min_threshold_k) == 0)
        elif method == 'k_bar':
            clr_condition = k_bar >= clr_threshold_condition
        else:
            raise ValueError(f"Invalid method {method} was specified")

        partly_cloudy_condition = k_bar >= partly_threshold_k

        if clr_condition:
            clr_date_list.append(date)
        elif partly_cloudy_condition:
            partly_cloudy_date_list.append(date)
        else:
            cloudy_date_list.append(date)
        
        average_k_list.append(k_bar)

    return clr_date_list, partly_cloudy_date_list, cloudy_date_list, average_k_list
