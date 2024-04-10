import pandas as pd
import numpy as np
import math 

from sklearn.linear_model import LinearRegression, HuberRegressor
from .clearsky_calculator import ClearskyCalculator
import matplotlib.pyplot as plt


class SkyConditionSplitter():
    """

    This class is designed to split any measurement data into three sky conditions: 
    clear-sky, partly-cloudy, and cloudy conditions. This is determined by the Rate Of Change (ROC), 
    which is defined by considering the differences in consecutive data point values. 
    For a day to be classified as clear-sky, it must have exactly one concave point, 
    and all data points on that day must be higher than the threshold_k. 
    For partly-cloudy conditions, if the day is not classified as clear-sky, 
    it will be considered partly-cloudy if the average k value for that day exceeds the splitting_k threshold. 
    Otherwise, it will be defined as a cloudy day.

    Parameters :

    measurement_df : pd.DataFrame
        The DataFrame should contain at least two fields: I (measured irradiance) and k (clear-sky index), calculated from measured irradiance divided by clear-sky irradiance.

    threshold_k : float
        The threshold that all data points in clear-sky day must pass.

    splitting_k : float
        The value that average clear-sky index of partly-cloudy data must pass.

    train_fac : float
        The proportion of training data after splitting.

    random_date : float
        The random seed of splitting.

    """

    
    def __init__(self, measurement_df, train_frac=0.6, threshold_k=0.7, splitting_k=0.6,random_state=42):
        """

        Initialize the instance by passing measurement data and other threshold values

        """
        self.measurement_df = measurement_df
        self.train_frac = train_frac
        self.random_state = random_state
        self.threshold_k = threshold_k
        self.splitting_k = splitting_k

        # Create dictionary attributes to store the number of samples for each condition in train and test sets
        self.train_samples_by_condition = {}
        self.test_samples_by_condition = {}

        # Call the methods during initialization to set 
        self.cluster_sky_condition()
        self.sample_train_test()

    @staticmethod
    def split_day_by_ROC(df_site, threshold_k=0.7, splitting_k=0.6):
        all_date_list = list(set(df_site.index.date))
        clr_date_list = []
        partly_cloudy_date_list = []
        cloudy_date_list = []

        for date in all_date_list:
            date_df = df_site[df_site.index.date == date].copy()  # Create a copy to avoid SettingWithCopyWarnings
            date_df['is_increasing'] = (date_df['I'].diff() > 0).astype(int)

            date_df['is_concave_point'] = date_df['is_increasing'].diff()
            date_df['is_concave_point'] = date_df['is_concave_point'] == -1

            clr_condition = (np.sum(date_df['is_concave_point']) == 1) and (np.sum(date_df['k'] <= threshold_k) == 0)
            partly_cloudy_condition = (date_df['k'].mean()) >= splitting_k

            if clr_condition:
                clr_date_list.append(date)
            elif partly_cloudy_condition:
                partly_cloudy_date_list.append(date)
            else:
                cloudy_date_list.append(date)

        return clr_date_list, partly_cloudy_date_list, cloudy_date_list

    
    def cluster_sky_condition(self, threshold_k =0.7, splitting_k=0.6):
        cluster_info_df = pd.DataFrame()

        for site_no in range(1, 57):
            site_name = f"site_{('00' + str(site_no))[-3:]}"
            df_site = self.measurement_df[self.measurement_df['site_name'] == site_name].copy()

            clr_date_list, partly_cloudy_date_list, cloudy_date_list = self.split_day_by_roc(df_site, threshold_k=threshold_k, spliting_k=splitting_k)

            cluster_df_site = pd.DataFrame({'Date': clr_date_list + partly_cloudy_date_list + cloudy_date_list,
                                            'site_name': [site_name] * len(clr_date_list + partly_cloudy_date_list + cloudy_date_list),
                                            'site_date': [f"{date} {site_name}" for date in clr_date_list + partly_cloudy_date_list + cloudy_date_list],
                                            'condition': ['clr']*len(clr_date_list) + ['partly_cloudy']*len(partly_cloudy_date_list) + ['cloudy']*len(cloudy_date_list)})
            cluster_info_df = pd.concat([cluster_info_df, cluster_df_site], ignore_index=True)
            print(f"Finished clustering sky-condition of site {site_name}")

        self.cluster_info_df = cluster_info_df

    def sample_train_test(self):
        train_df = pd.DataFrame()
        test_df = pd.DataFrame()

        for condition in ['clr', 'partly_cloudy', 'cloudy']:
            condition_df = self.cluster_info_df[self.cluster_info_df['condition'] == condition]
            condition_samples = condition_df['site_date'].tolist()

            # Filter the original DataFrame based on the selected samples for the current condition
            condition_data = self.measurement_df[self.measurement_df['site_date'].isin(condition_samples)]

            # Split the filtered data into train and test sets
            condition_train_df = condition_data.sample(frac=self.train_frac, random_state=self.random_state)
            condition_test_df = condition_data.drop(condition_train_df.index)

            train_df = pd.concat([train_df, condition_train_df], ignore_index=True)
            test_df = pd.concat([test_df, condition_test_df], ignore_index=True)

            # Store the number of samples for each condition in train and test sets
            self.train_samples_by_condition[condition] = len(condition_train_df)
            self.test_samples_by_condition[condition] = len(condition_test_df)
        self.train_df = train_df
        self.test_df = test_df

    @staticmethod
    def cal_new_TL(clr_data_df, alt, huber_epsilon=1):
        _df = clr_data_df.copy()
        I0 = 1366.1
        h = alt
        f1 = math.exp(-h/8000);
        f2 = math.exp(-h/1250);
        a1 = (h*5.09e-5) + 0.868;
        a2 = (h*3.92e-5) + 0.0387;
        dependent_var =  ((f1 - f2)*a2*_df['AM'] + (np.log(_df['I'] / (a1 * I0 * _df['cos_zna'])))).to_numpy().reshape(-1, 1)
        independent_var = (- a2 * f2 * _df['AM']).to_numpy().reshape(-1, 1)
        
        linear_model = LinearRegression()
        huber_loss_model = HuberRegressor(epsilon=huber_epsilon)

        linear_model.fit(independent_var, dependent_var)
        huber_loss_model.fit(independent_var, dependent_var)

        TL_linear = linear_model.coef_[0]
        TL_huber = huber_loss_model.coef_[0]

        return TL_linear[0], TL_huber

    @staticmethod
    def update_new_clr(df_site, lat, long, new_TL):

        _df = df_site.copy()

        _start_date = _df.index.date[0].strftime('%Y-%m-%d')
        _end_date = _df.index.date[-1].strftime('%Y-%m-%d')

        _site_obj = ClearskyCalculator(lat=lat, long=long, defined_TL=new_TL)

        _iclr_df = _site_obj.cal_clearsky_irradiance(start_date=_start_date, end_date=_end_date)

        _merged_df = pd.merge(_df, _iclr_df[['Iclr']], left_index=True, right_index=True, how='left', suffixes=('_df', '_iclr_df'))

        _df['Iclr'] = _merged_df['Iclr_iclr_df'].values
        _df.loc[:, 'k'] = _df.loc[:, 'I'] / _df.loc[:, 'Iclr']

        return _df
    
    @staticmethod
    def plot_updated_irradiance(old_df_site, clr_date_list, lat, long, alt, ncols=3):
        nrows = (len(clr_date_list) // ncols) +1

        fig, ax = plt.subplots(nrows,ncols, figsize=(10*ncols, 10*nrows))

        # Reshape to access by 1D-array(row) of ax
        ax = ax.reshape(1 ,-1)[0]

        _clr_data_df = old_df_site[np.isin(old_df_site.index.date, clr_date_list)]

        _, new_TL = SkyConditionSplitter.cal_new_TL(_clr_data_df, alt, huber_epsilon=1)
        _updated_df_site = SkyConditionSplitter.update_new_clr(old_df_site, lat=lat, long=long, new_TL=new_TL)


        for k in range(len(ax)):
            plt.subplot(nrows, ncols, k+1)

            if k >= len(clr_date_list):
                break
            plot_date = clr_date_list[k]
            old_date_df = old_df_site[old_df_site.index.date==plot_date]
            new_date_df = _updated_df_site[_updated_df_site.index.date==plot_date]

            line1 = plt.plot(old_date_df['I'], '-x', color='black')
            line2 = plt.plot(old_date_df['Iclr'], '-^', color='red')
            line3 = plt.plot(new_date_df['Iclr'], '-o', color='green')

            if k%ncols ==0:
                plt.ylabel('Irradiance [W/m2]', fontsize=20)
            plt.title(plot_date,fontsize=20)
            plt.xticks(rotation = 'vertical')
            plt.xlabel('Time', fontsize=20)
            plt.grid(True)
            plt.subplots_adjust(hspace=0.9)

        fig.legend(labels=['I','old_Iclr', 'updated_Iclr'],fontsize=20)

        return fig



        





        
