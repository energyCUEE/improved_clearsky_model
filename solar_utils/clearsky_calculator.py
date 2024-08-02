import pandas as pd
import numpy as np
import os
from typing import Optional, Tuple, Union

class ClearskyCalculator:
    """
    This class provides methods to calculate clear-sky irradiance variables based on the Ineichen model. 
    Additionally, it offers the flexibility to adjust the Linke turbidity coefficient (TL) according to specific locations,
    also with informed by the monthly TL values derived from the research presented in the Janjai paper,
    and also the defined monthly TL that users can specify through the constructor.

    Janjai, S., Kumharn, W., & Laksanaboonsong, J. (2003). Determination of Angstrom's turbidity coefficient over Thailand. 
    Renewable Energy, 28(11), 1685-1700.
    """

    # Locate the metadata directory
    _module_dir = os.path.dirname(__file__)
    _matadata_dir = os.path.join(_module_dir, "metadata")
    
    # Read all relevant metadata 
    _gridtiff_metadata_df = pd.read_csv(os.path.join(_matadata_dir, 'plant_metadata.csv'))
    _lat_metadata_arr = (_gridtiff_metadata_df['latt'] * np.pi / 180).to_numpy().reshape(-1, 1)
    _long_metadata_arr = (_gridtiff_metadata_df['long'] * np.pi / 180).to_numpy().reshape(-1, 1)
    _alt_metadata_arr = _gridtiff_metadata_df['alt'].to_numpy()
    _solarmap_TL_metadata_arr = _gridtiff_metadata_df['old_TL'].to_numpy()
    
    _JanJai_metadata_df = pd.read_csv(os.path.join(_matadata_dir, 'JanJaiAngstromTurbidity.csv'))
    _JanJai_lat_metadata_arr = (_JanJai_metadata_df['Latitude'] * np.pi / 180).to_numpy().reshape(-1, 1)
    _JanJai_long_metadata_arr = (_JanJai_metadata_df['Longitude'] * np.pi / 180).to_numpy().reshape(-1, 1)
    
    def __init__(self, lat: float, long: float, defined_TL: Optional[float] = None, 
                 defined_alt: Optional[float] = None, monthly_TL: Optional[np.ndarray] = None):
        """
        Initialize an instance of ClearskyCalculator with relevant attributes.

        Parameters:
        lat (float): Latitude of the location in degrees.
        long (float): Longitude of the location in degrees.
        defined_TL (Optional[float]): User-defined fixed Linke turbidity coefficient for all months. Default is None.
        defined_alt (Optional[float]): User-defined altitude of the location in meters. Default is None.
        monthly_TL (Optional[np.ndarray]): User-defined array of monthly Linke turbidity coefficients. Default is None.

        Attributes:
        lat (float): Latitude of the location in degrees.
        long (float): Longitude of the location in degrees.
        nearest_alt (float): Nearest altitude value from metadata based on the given location.
        nearest_estimated_TL (float): Nearest estimated TL value from metadata based on the given location of all time.
        JanJai_monthly_TL (np.ndarray): Monthly TL values derived from Janjai's research for the nearest location.
        defined_TL (Optional[float]): User-defined fixed TL value.
        defined_alt (Optional[float]): User-defined altitude value.
        monthly_TL (Optional[np.ndarray]): User-defined array of monthly TL values.
        """
        self.lat = lat
        self.long = long
        self.nearest_alt, self.nearest_estimated_TL = self.get_nearest_measurement_metadata()
        self.JanJai_monthly_TL = self.get_nearest_JanJai_metadata()
        self.defined_TL = defined_TL
        self.defined_alt = defined_alt
        self.monthly_TL = monthly_TL.ravel() if monthly_TL is not None else None
        
    def get_nearest_measurement_metadata(self) -> Tuple[float, float]:
        """
        Determine the nearest metadata by calculating the spherical distance based on latitude and longitude.

        Returns:
        Tuple[float, float]: Nearest altitude and estimated TL value.
        """
        user_rad_lat = self.lat * np.pi / 180
        user_rad_long = self.long * np.pi / 180
        dist_arr = np.arccos(np.sin(user_rad_lat) * np.sin(ClearskyCalculator._lat_metadata_arr) +
                             np.cos(user_rad_lat) * np.cos(ClearskyCalculator._lat_metadata_arr) * 
                             np.cos(abs(user_rad_long - ClearskyCalculator._long_metadata_arr)))
        nearest_index = np.argmin(dist_arr)
        nearest_alt = ClearskyCalculator._alt_metadata_arr[nearest_index]
        nearest_estimated_TL = ClearskyCalculator._solarmap_TL_metadata_arr[nearest_index]
        return nearest_alt, nearest_estimated_TL
        
    def get_nearest_JanJai_metadata(self) -> np.ndarray:
        """
        Retrieve the nearest Janjai metadata based on spherical distance calculation.

        Returns:
        np.ndarray: Monthly TL values derived from Janjai's research for the nearest location.
        """
        user_rad_lat = self.lat * np.pi / 180
        user_rad_long = self.long * np.pi / 180
        dist_arr = np.arccos(np.sin(user_rad_lat) * np.sin(ClearskyCalculator._JanJai_lat_metadata_arr) +
                             np.cos(user_rad_lat) * np.cos(ClearskyCalculator._JanJai_lat_metadata_arr) * 
                             np.cos(abs(user_rad_long - ClearskyCalculator._JanJai_long_metadata_arr)))
        nearest_index = np.argmin(dist_arr)
        chosen_beta_df = ClearskyCalculator._JanJai_metadata_df.loc[nearest_index, 
                                                                    ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 
                                                                     'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']]
        
        AngstromTurbidity = chosen_beta_df.to_numpy()
        
        # Precipitable Water (cm) downloaded from AERONET (averaged over 2019-2022) at BKK lat/long
        PW = np.array([3.0429, 3.3528, 4.0548, 4.5590, 5.2208, 5.5694, 4.7205, 5.2554, 5.1380, 4.7121, 4.2269, 3.2973])
        
        # Remund Formula of TL for each month
        month_TL = (1.8494 + 0.24245 * PW - 0.0203 * PW**2) + (15.427 + 0.3153 * PW - 0.0254 * PW**2) * AngstromTurbidity
        return month_TL
        
    def cal_cos_zna(self, start_date: str, end_date: str, freq: str = '15min') -> pd.DataFrame:
        """
        Calculate cosine of solar zenith angle.

        Parameters:
        start_date (str): Start date in 'YYYY-MM-DD' format.
        end_date (str): End date in 'YYYY-MM-DD' format.
        freq (str): Frequency of time intervals. Default is '15min'.

        Returns:
        pd.DataFrame: DataFrame containing datetime and cosine of solar zenith angle.
        """
        phi = self.lat
        date_range = pd.date_range(start=f"{start_date} 00:00", end=f"{end_date} 23:59", freq=freq)
        doy_arr = date_range.dayofyear.to_numpy()
        hour_arr = date_range.hour.to_numpy()
        minute_arr = date_range.minute.to_numpy()

        time_arr = hour_arr + minute_arr / 60
        B = (360 / 365) * (doy_arr - 81)
        LSTM = 15 * 7
        EoT = 9.87 * np.sin(2 * B * np.pi / 180) - 7.53 * np.cos(B * np.pi / 180) - 1.5 * np.sin(B * np.pi / 180)
        TC = 4 * (self.long - LSTM) + EoT
        LST = time_arr + TC / 60
        HRA = 15 * abs(LST - 12)
        delta = -23.45 * np.cos(2 * np.pi * (doy_arr + 10) / 365)
        cos_zna = np.sin(delta * np.pi / 180) * np.sin(phi * np.pi / 180) + np.cos(delta * np.pi / 180) * np.cos(phi * np.pi / 180) * np.cos(HRA * np.pi / 180)
        
        cos_zna_df = pd.DataFrame({'Datetime': date_range, 'cos_zna': cos_zna})
        cos_zna_df.set_index('Datetime', inplace=True)
        
        return cos_zna_df
    
    def cal_dairmass(self, start_date: str, end_date: str, freq: str = '15min', 
                     cos_zna_df: Optional[pd.DataFrame] = None) -> pd.DataFrame:
        """
        Calculate airmass.

        Parameters:
        start_date (str): Start date in 'YYYY-MM-DD' format.
        end_date (str): End date in 'YYYY-MM-DD' format.
        freq (str): Frequency of time intervals. Default is '15min'.
        cos_zna_df (Optional[pd.DataFrame]): DataFrame containing cosine of solar zenith angle. Default is None.

        Returns:
        pd.DataFrame: DataFrame containing datetime and airmass.
        """
        if cos_zna_df is None:
            cos_zna_df = self.cal_cos_zna(start_date, end_date, freq)
        cos_zna = cos_zna_df['cos_zna'].to_numpy()
        cos_zna = np.where(cos_zna < 0, 0, cos_zna)
        AM_arr = 1 / (cos_zna + 0.50572 * (96.07995 - np.degrees(np.arccos(cos_zna))) ** -1.6364)  # Kasten and Young
        AM_df = pd.DataFrame({'Datetime': cos_zna_df.index, 'AM': AM_arr})
        AM_df.set_index('Datetime', inplace=True)
        return AM_df
    
    def cal_clearsky_irradiance(self, start_date: str, end_date: str, freq: str = '15min', 
                                choice: str = 'monthly_estimate', cos_zna_df: Optional[pd.DataFrame] = None, 
                                AM_df: Optional[pd.DataFrame] = None) -> pd.DataFrame:
        """
        Calculate clear-sky irradiance using a specific TL by choice.

        Parameters:
        start_date (str): Start date in 'YYYY-MM-DD' format.
        end_date (str): End date in 'YYYY-MM-DD' format.
        freq (str): Frequency of time intervals. Default is '15min'.
        choice (str): Choice of TL value source ('monthly_estimate', 'estimate', 'JanJai'). Default is 'monthly_estimate'.
        cos_zna_df (Optional[pd.DataFrame]): DataFrame containing cosine of solar zenith angle. Default is None.
        AM_df (Optional[pd.DataFrame]): DataFrame containing airmass. Default is None.

        Returns:
        pd.DataFrame: DataFrame containing datetime and clear-sky irradiance.
        """
        if cos_zna_df is None:
            cos_zna_df = self.cal_cos_zna(start_date, end_date, freq)
        if AM_df is None:
            AM_df = self.cal_dairmass(start_date, end_date, freq, cos_zna_df)
        
        cos_zna = cos_zna_df['cos_zna'].to_numpy()
        AM_arr = AM_df['AM'].to_numpy()
        
        date_range = pd.date_range(start=f"{start_date} 00:00", end=f"{end_date} 23:59", freq=freq)
        
        TL_arr = np.zeros(len(date_range))

        if self.defined_TL is not None:
            TL_arr[:] = self.defined_TL
        elif choice == 'estimate':
            TL_arr[:] = self.nearest_estimated_TL
        elif choice == 'JanJai':
            month_arr = date_range.month.to_numpy()
            for month_index in range(12):
                TL_arr = TL_arr + np.equal(month_arr, month_index + 1) * self.JanJai_monthly_TL[month_index]
        elif choice == 'monthly_estimate':
            month_arr = date_range.month.to_numpy()
            for month_index in range(12):
                TL_arr = TL_arr + np.equal(month_arr, month_index + 1) * self.monthly_TL[month_index]
        else:
            raise ValueError(f"Invalid choice {choice} is specified. Please specify a valid choice.")

        if self.defined_TL is not None:
            alt = self.defined_alt
        else:
            alt = self.nearest_alt

        f1 = np.e ** (-alt / 8000)
        f2 = np.e ** (-alt / 1250)
        a1 = (alt * 5.09e-5) + 0.868
        a2 = (alt * 3.92e-5) + 0.0387
        Isc = 1366.1
        iclr_arr = a1 * Isc * cos_zna * np.e ** (-a2 * AM_arr * (f1 + f2 * (TL_arr - 1)))
        iclr_arr = np.clip(iclr_arr, a_min=0, a_max=None)

        iclr_df = pd.DataFrame({'Datetime': date_range, 'Iclr': iclr_arr})
        iclr_df.set_index('Datetime', inplace=True)
        return iclr_df
        
    def get_solar_info(self, start_date: str, end_date: str, freq: str = '15min', choice: str = 'monthly_estimate') -> pd.DataFrame:
        """
        Retrieve all calculated solar information in DataFrame format.

        Parameters:
        start_date (str): Start date in 'YYYY-MM-DD' format.
        end_date (str): End date in 'YYYY-MM-DD' format.
        freq (str): Frequency of time intervals. Default is '15min'.
        choice (str): Choice of TL value source ('monthly_estimate', 'estimate', 'JanJai'). Default is 'monthly_estimate'.

        Returns:
        pd.DataFrame: DataFrame containing datetime, cosine of solar zenith angle, airmass, and clear-sky irradiance.
        """
        cos_zna_df = self.cal_cos_zna(start_date, end_date, freq)
        AM_df = self.cal_dairmass(start_date, end_date, freq, cos_zna_df)
        Iclr_df = self.cal_clearsky_irradiance(start_date, end_date, freq, choice=choice, cos_zna_df=cos_zna_df, AM_df=AM_df)
        
        solar_info_df = cos_zna_df.join(AM_df, how='inner').join(Iclr_df, how='inner')
        
        return solar_info_df
