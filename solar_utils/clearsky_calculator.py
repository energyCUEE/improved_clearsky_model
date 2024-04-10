import pandas as pd
import numpy as np
import os


class ClearskyCalculator:

    """
    This class provides methods to calculate clear-sky irradiance variables based on the Ineichen model. 
    Additionally, it offers the flexibility to adjust the Linke turbidity coefficient (TL) according to specific locations.
    also with informed by the TL values derived from 
    the research presented in the Janjai paper.

    Janjai, S., Kumharn, W., & Laksanaboonsong, J. (2003). Determination of Angstrom's turbidity coefficient over Thailand. 
    Renewable Energy, 28(11), 1685-1700.

    
    """

    # Locate the metadata directory

    _module_dir = os.path.dirname(__file__)
    _matadata_dir = os.path.join(_module_dir, "metadata")
    
    # Read all relevant metadata 

    _gridtiff_metadata_df = pd.read_csv(os.path.join(_matadata_dir,'plant_metadata.csv'))
    _lat_metadata_arr = (_gridtiff_metadata_df['lat']*np.pi/180).to_numpy().reshape(-1,1)
    _long_metadata_arr = (_gridtiff_metadata_df['long']*np.pi/180).to_numpy().reshape(-1,1)
    
    _alt_metadata_arr = _gridtiff_metadata_df['alt'].to_numpy()
    _solarmap_TL_metadata_arr = _gridtiff_metadata_df['TL'].to_numpy()
    
    _JanJai_metadata_df = pd.read_csv(os.path.join(_matadata_dir,'JanJaiAngstromTurbidity.csv'))
    _JanJai_lat_metadata_arr = (_JanJai_metadata_df['Latitude']*np.pi/180).to_numpy().reshape(-1,1)
    _JanJai_long_metadata_arr = (_JanJai_metadata_df['Longitude']*np.pi/180).to_numpy().reshape(-1,1)
    
    

    def __init__(self, lat, long, defined_TL=None):
        """
        Initialize an instance with latitude and longitude, You can use your own TL by specify fix TL in defined_TL
        
        """
        self.lat = lat
        self.long = long
        self.nearest_alt, self.nearest_estimated_TL = self.get_nearest_measurement_metadata()
        self.month_TL = self.get_nearest_JanJai_metadata()
        self.defined_TL = defined_TL
        
    # Determine the nearest metadata by calculating the spherical distance based on latitude and longitude
        
    def get_nearest_measurement_metadata(self):
        user_rad_lat = self.lat*np.pi/180
        user_rad_long = self.long*np.pi/180
        dist_arr = np.arccos(np.sin(user_rad_lat)*np.sin(ClearskyCalculator._lat_metadata_arr)+
                     np.cos(user_rad_lat)*np.cos(ClearskyCalculator._lat_metadata_arr )*np.cos(abs(user_rad_long-ClearskyCalculator._long_metadata_arr)))
        nearest_index = np.argmin(dist_arr)
        nearest_alt = ClearskyCalculator._alt_metadata_arr[nearest_index]
        nearest_estimated_TL = ClearskyCalculator._solarmap_TL_metadata_arr[nearest_index]
        return nearest_alt, nearest_estimated_TL
        
    def get_nearest_JanJai_metadata(self):
        user_rad_lat = self.lat*np.pi/180
        user_rad_long = self.long*np.pi/180
        dist_arr = np.arccos(np.sin(user_rad_lat)*np.sin(ClearskyCalculator._JanJai_lat_metadata_arr)+
                     np.cos(user_rad_lat)*np.cos(ClearskyCalculator._JanJai_lat_metadata_arr )*np.cos(abs(user_rad_long-ClearskyCalculator._JanJai_long_metadata_arr)))
        nearest_index = np.argmin(dist_arr)
        chosen_beta_df = ClearskyCalculator._JanJai_metadata_df.loc[nearest_index, ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']]
        
        AngstromTurbidity = chosen_beta_df.to_numpy()
        
        # Precitable Water (cm) were downloaded from AERONET (averaged over 2019-2022) at BKK lat/long
        PW = np.array( [3.0429, 3.3528, 4.0548, 4.5590, 5.2208, 5.5694, 4.7205, 5.2554, 5.1380, 4.7121, 4.2269, 3.2973 ] )
        # Remund Formula of TL for each month -- returned as 12-entry array

        month_TL = (1.8494+0.24245 * PW - 0.0203 * PW**2) + (15.427+0.3153 * PW - 0.0254 * PW ** 2) * AngstromTurbidity
        return month_TL
        
        
    # Calculate solar zenith angle in degree units
    
    def cal_dzna(self, start_date, end_date, freq='15min'):
        phi = self.lat
        date_range = pd.date_range(start=start_date, end=end_date, freq=freq)
        doy_arr = date_range.dayofyear.to_numpy()
        hour_arr = date_range.hour.to_numpy()
        minute_arr = date_range.minute.to_numpy()

        time_arr = hour_arr + minute_arr/60
        B = (360/365)*(doy_arr-81)
        LSTM = 15*7
        EoT = 9.87*np.sin(2*B*np.pi/180)-7.53*np.cos(B*np.pi/180)-1.5*np.sin(B*np.pi/180)
        TC = 4*(self.long-LSTM)+EoT
        LST = time_arr + TC/60
        HRA = 15*abs(LST-12)
        delta= -23.45*np.cos(2*np.pi*(doy_arr+10)/365)
        zna = np.sin(delta*np.pi/180)*np.sin(phi*np.pi/180) + np.cos(delta*np.pi/180)*np.cos(phi*np.pi/180)*np.cos(HRA*np.pi/180)
        zna_arr = np.arccos(zna)*180/np.pi
        
        zna_df = pd.DataFrame({'Datetime':date_range, 'zna':zna_arr})
        zna_df.set_index('Datetime', inplace=True)
        
        return zna_df
    
    # Calculate airmass
    def cal_dairmass(self, start_date, end_date, freq='15min'):
        phi = self.lat
        date_range = pd.date_range(start=start_date, end=end_date, freq=freq)
        doy_arr = date_range.dayofyear.to_numpy()
        hour_arr = date_range.hour.to_numpy()
        minute_arr = date_range.minute.to_numpy()
        time_arr = hour_arr + minute_arr/60
        B = (360/365)*(doy_arr - 81)
        LSTM = 15*7
        EoT = 9.87*np.sin(2*B*np.pi/180)-7.53*np.cos(B*np.pi/180)-1.5*np.sin(B*np.pi/180)
        TC = 4*(self.long-LSTM)+EoT
        LST = time_arr + TC/60
        HRA = 15*abs(LST-12)
        delta= -23.45*np.cos(2*np.pi*(doy_arr+10)/365)
        x = np.sin(delta*np.pi/180)*np.sin(phi*np.pi/180) + np.cos(delta*np.pi/180)*np.cos(phi*np.pi/180)*np.cos(HRA*np.pi/180)
        x = np.where(x<0, 0, x)
        AM_arr = 1 / (x + 0.50572 * (96.07995 - np.degrees(np.arccos(x))) ** -1.6364) # Kasten and Young
        AM_df = pd.DataFrame({'Datetime':date_range, 'AM':AM_arr})
        AM_df.set_index('Datetime', inplace=True)
        return AM_df
    
        
    # Calculate clear-sky irradiance using specific TL by choice

    def cal_clearsky_irradiance(self, start_date, end_date, freq='15min', choice='estimate'):
        phi = self.lat
        date_range = pd.date_range(start=start_date, end=end_date, freq=freq)
        doy_arr = date_range.dayofyear.to_numpy()
        hour_arr = date_range.hour.to_numpy()
        minute_arr = date_range.minute.to_numpy()
        time_arr = hour_arr + minute_arr/60
        B = (360/365)*(doy_arr - 81)
        LSTM = 15*7
        EoT = 9.87*np.sin(2*B*np.pi/180)-7.53*np.cos(B*np.pi/180)-1.5*np.sin(B*np.pi/180)
        TC = 4*(self.long-LSTM)+EoT
        LST = time_arr + TC/60
        HRA = 15*abs(LST-12)
        delta= -23.45*np.cos(2*np.pi*(doy_arr+10)/365)
        x = np.sin(delta*np.pi/180)*np.sin(phi*np.pi/180) + np.cos(delta*np.pi/180)*np.cos(phi*np.pi/180)*np.cos(HRA*np.pi/180)
        x = np.where(x<0, 0, x)
        AM_arr = 1 / (x + 0.50572 * (96.07995 - np.degrees(np.arccos(x))) ** -1.6364)
        
        TL_arr = np.zeros(len(date_range))

        if self.defined_TL is not None:
            TL_arr[:] = self.defined_TL
        elif choice == 'estimate':
            TL_arr[:] = self.nearest_estimated_TL
        elif choice == 'JanJai':
            month_arr = date_range.month.to_numpy()
            for month_index in range(12):
                TL_arr = TL_arr + np.equal(month_arr, month_index+1)*self.month_TL[month_index]
        else:
            raise ValueError(f"Invalid choice {choice} is specified. Please specify a valid choice.")
                
        f1 = np.e ** (-self.nearest_alt / 8000)
        f2 = np.e ** (-self.nearest_alt / 1250)
        a1 = (self.nearest_alt * 5.09e-5) + 0.868
        a2 = (self.nearest_alt * 3.92e-5) + 0.0387
        Isc = 1366.1
        iclr_arr = a1 * Isc * x * np.e ** (-a2 * AM_arr * (f1 + f2 * (TL_arr - 1)))
        iclr_arr = np.clip(iclr_arr, a_min=0, a_max=None)
        
        iclr_df = pd.DataFrame({'Datetime':date_range, 'Iclr':iclr_arr})
        iclr_df.set_index('Datetime', inplace=True)
        return iclr_df
        
    # Retrieve all calculated information in DataFrame format

    def get_solar_info(self, start_date, end_date, freq='15min', choice='estimate'):
        zna_df = self.cal_dzna(start_date, end_date, freq)
        AM_df = self.cal_dairmass(start_date, end_date, freq)
        Iclr_df = self.cal_clearsky_irradiance(start_date, end_date, freq, choice)
        
        solar_info_df = zna_df.join(AM_df, how='inner').join(Iclr_df, how='inner')
        solar_info_df['cos_zna'] = np.cos(solar_info_df['zna']*np.pi/180)
        
        return solar_info_df
    
   