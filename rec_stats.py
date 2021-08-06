from rec import TimebinInterval
from versioning import *
from rec_db import *
from sklearn.preprocessing import StandardScaler

class PopulationBehavioralTimeseries():

    base_columns = ['Unit', 'Event', 'Condition', 'Instance', 'Timeseries']

    def __init__(self, condition_labels=[], timebin_interval=TimebinInterval(0, 0, 0, 0)):
        self.df = pd.DataFrame(columns=self.base_columns)
        self.condition_labels = condition_labels
        self.timebin_interval = timebin_interval

    def init_with_df(self, df):
        pbt = PopulationBehavioralTimeseries(condition_labels=self.condition_labels, timebin_interval=self.timebin_interval)
        pbt.df = df
        return pbt

    def get_unit_inds(self):
        return list(self.df['Unit'].unique())

    def get_condition_levels_1d(self):
        return list(self.df['Condition'].unique())

    def get_condition_levels_nd(self):
        cond_1d_list = self.get_condition_levels_1d()
        return [tuple(np.unique(t)) for t in zip(*cond_1d_list)]

    def get_condition_slice(self, condition_list):
        return self.df.loc[self.df['Condition'].isin(condition_list)]

    def get_unit_slice(self, unit_list):
        return self.df.loc[self.df['Unit'].isin(unit_list)]

    def add_data_row(self, unit_ind, event_ind, condition, instance, timeseries, index=None):
        row_data_df = pd.DataFrame.from_dict({index: [unit_ind, event_ind, condition, instance, timeseries]},
                                             columns=self.base_columns,
                                             orient='index')
        ignore_index = index is None
        self.df.append(row_data_df, ignore_index=ignore_index)

    def add_data_rows_from_list(self, data_list):
        row_data_df = pd.DataFrame(data_list, columns=self.base_columns)
        self.df = self.df.append(row_data_df, ignore_index=False)

    def func_instances(self, group_columns, func):
        grouper = self.df.groupby(group_columns)
        grouper_index = grouper.groups.keys()
        timeseries = grouper['Timeseries'].aggregate(func)
        events = grouper['Event'].aggregate(list)
        df_average = pd.DataFrame({'Timeseries': timeseries,
                                   'Event': events,
                                   'Instance': [0 for _ in range(len(grouper_index))]},
                                  index=grouper_index)
        df_average.index.set_names(group_columns, inplace=True)
        df_average.reset_index(drop=False, inplace=True)
        return self.init_with_df(df_average)

    def average_instances(self, group_columns):
        return self.func_instances(group_columns, lambda x: list(np.mean(np.array(list(x)), axis=0)))

    def std_instances(self, group_columns):
        return self.func_instances(group_columns, lambda x: list(np.std(np.array(list(x)), axis=0)))

    def se_instances(self, group_columns):
        return self.func_instances(group_columns, lambda x: list(np.std(np.array(list(x)), axis=0) / np.sqrt(len(x))))

    def crop_timeseries(self, t_start, t_end):
        crop_start_ind = self.timebin_interval.split_to_bins_offset().index(t_start)
        crop_end_ind = self.timebin_interval.split_to_bins_offset().index(t_end)
        crop_slice = slice(crop_start_ind, crop_end_ind + 1)
        self.df['Timeseries'] = self.df['Timeseries'].apply(lambda x: x[crop_slice])
        self.timebin_interval = self.timebin_interval.sub_interval(t_start - self.timebin_interval.timebin, t_end)

    def smooth_timeseries(self, smoother):
        return self.df['Timeseries'].apply(smoother.smoothen)

    def smooth_df(self, smoother):
        smooth_timeseries = self.smooth_timeseries(smoother)
        df_smooth = self.df.copy()
        df_smooth['Timeseries'] = smooth_timeseries
        return df_smooth

    def unfold_conditions_df(self):
        return pd.DataFrame(self.df['Condition'].tolist(), index=self.df.index, columns=self.condition_labels)

    def to_PCA_array(self):
        df_base = self.df
        grouper = df_base.groupby(['Unit', 'Condition', 'Instance'])
        num_units = df_base['Unit'].nunique()
        num_conditions = df_base['Condition'].nunique()
        num_instances = df_base.groupby(['Unit', 'Condition']).size().unique()[0]
        num_timebins = self.timebin_interval.num_of_bins()
        # order indices by groupping order unit, condition, instance
        index_list = [index for _, group in grouper for index in group.index]
        # create pre_pca unit, condition, instance, timebin array
        X_preshape = (num_units, num_conditions, num_instances, num_timebins)
        X_pre = np.array(df_base.iloc[index_list]['Timeseries'].tolist()).reshape(X_preshape)
        # convert to pca array format
        X_shape = (num_conditions * num_instances * num_timebins, num_units)
        X = X_pre.transpose(1, 2, 3, 0).reshape(X_shape)
        # unit, condition, instance information for reconstruction
        records = df_base.iloc[index_list][['Unit', 'Condition', 'Instance']].to_records()
        return X, records

    def to_PCA_scale_array(self):
        X, records = self.to_PCA_array()
        X_scale = StandardScaler().fit_transform(X.transpose()).transpose()
        return (X, X_scale), records

    def to_dPCA_mean_array(self):
        # init params
        pbt_mean = self.average_instances(['Unit', 'Condition'])
        df_base = pd.concat([pbt_mean.df, pbt_mean.unfold_conditions_df()], axis=1)
        grouper = df_base.groupby(['Unit'] + pbt_mean.condition_labels + ['Instance'])
        num_units = df_base['Unit'].nunique()
        num_conditions = list(df_base[pbt_mean.condition_labels].nunique())
        num_instances = df_base.groupby(['Unit', 'Condition']).size().unique()[0]
        num_timebins = pbt_mean.timebin_interval.num_of_bins()
        # order indices by groupping order unit, condition, instance
        index_list = [index for _, group in grouper for index in group.index]
        # create pre_pca unit, condition, instance, timebin array
        X_preshape = tuple([num_units] + num_conditions + [num_instances, num_timebins])
        X_pre = np.array(df_base.iloc[index_list]['Timeseries'].tolist()).reshape(X_preshape)
        # convert to pca array format
        X = X_pre.squeeze()
        # unit, condition, instance information for reconstruction
        records = df_base.iloc[index_list][['Unit', 'Condition', 'Instance']].to_records()
        return X, records

    def to_dPCA_mean_demean_array(self):
        X_mean, records = self.to_dPCA_mean_array()
        mean_shape = X_mean.shape
        X__mean_2d = X_mean.reshape((mean_shape[0], -1))
        X_mean_demean = StandardScaler(with_std=False).fit_transform(X__mean_2d.transpose()).transpose().reshape(mean_shape)
        return (X_mean, X_mean_demean), records

    def to_dPCA_trial_array(self):
        # init params
        df_base = pd.concat([self.df, self.unfold_conditions_df()], axis=1)
        grouper = df_base.groupby(['Unit'] + self.condition_labels + ['Instance'])
        num_units = df_base['Unit'].nunique()
        num_conditions = list(df_base[self.condition_labels].nunique())
        num_instances = df_base.groupby(['Unit', 'Condition']).size().unique()[0]
        num_timebins = self.timebin_interval.num_of_bins()
        # order indices by groupping order unit, condition, instance
        index_list = [index for _, group in grouper for index in group.index]
        # create pre_pca unit, condition, instance, timebin array
        X_preshape = tuple([num_units] + num_conditions + [num_instances, num_timebins])
        X_pre = np.array(df_base.iloc[index_list]['Timeseries'].tolist()).reshape(X_preshape)
        # convert to pca array format
        X_transpose_order = tuple([len(num_conditions) + 1, 0] + list(range(1, len(num_conditions) + 1)) + [len(num_conditions) + 2])
        X = X_pre.transpose(X_transpose_order)
        # unit, condition, instance information for reconstruction
        records = df_base.iloc[index_list][['Unit', 'Condition', 'Instance']].to_records()
        return X, records

    def to_pseudosession_firing_rate(self, pseudoarray_inds, t_ind):

        df_ordered_timeseries = self.df.loc[pseudoarray_inds.flat]['Timeseries']
        df_t_fr = df_ordered_timeseries.apply(lambda timeseries: timeseries[t_ind])
        X = df_t_fr.to_numpy().reshape(pseudoarray_inds.shape)
        return X

class FactorBehavioralTimeseries(PopulationBehavioralTimeseries):

    base_columns = ['Factor', 'Condition', 'Instance', 'Timeseries']

    def __init__(self, df, condition_labels, timebin_interval):
        self.df = df
        self.condition_labels = condition_labels
        self.timebin_interval = timebin_interval

    def init_with_df(self, df):
        pbt = PopulationBehavioralTimeseries(condition_labels=self.condition_labels,
                                             timebin_interval=self.timebin_interval)
        pbt.df = df
        return pbt

    def unfold_conditions_df(self):
        return pd.DataFrame(self.df['Condition'].tolist(), index=self.df.index, columns=self.condition_labels)

    def average_instances(self, group_columns):
        grouper = self.df.groupby(group_columns)
        grouper_index = grouper.groups.keys()
        timeseries = grouper['Timeseries'].aggregate(lambda x: list(np.mean(np.array(list(x)), axis=0)))
        df_average = pd.DataFrame({'Timeseries': timeseries,
                                   'Instance': [0 for _ in range(len(grouper_index))]},
                                  index=grouper_index)
        df_average.index.set_names(group_columns, inplace=True)
        df_average.reset_index(drop=False, inplace=True)
        return FactorBehavioralTimeseries(df_average, self.condition_labels, self.timebin_interval)

    def init_with_df(self): pass
    def get_unit_inds(self): pass
    def get_unit_slice(self): pass
    def add_data_row(self): pass
    def add_data_rows_from_list(self): pass
    def to_PCA_array(self): pass
    def to_dPCA_mean_array(self): pass
    def to_dPCA_trial_array(self): pass

    def fbt_df_from_PCA(X_factor, records, decomp_obj, timebin_interval):
        records_df = pd.DataFrame(records)
        num_factors = decomp_obj.n_components_
        num_conditions = records_df['Condition'].nunique()
        num_instances = records_df['Instance'].nunique()
        num_timebins = timebin_interval.num_of_bins()
        X_factor_base = X_factor.transpose().reshape(num_factors, num_conditions, num_instances, num_timebins)
        records_base = list(
            product(list(range(num_factors)), list(records_df['Condition'].unique()), list(range(num_instances))))
        records_base_df = pd.DataFrame(records_base, columns=['Factor', 'Condition', 'Instance'])
        records_base_df['Timeseries'] = list(X_factor_base.reshape(len(records_base_df), num_timebins))
        return records_base_df

    def fbt_df_from_dPCA(X_factor, records, decomp_obj, timebin_interval):
        records_df = pd.DataFrame(records)
        num_factors = decomp_obj.n_components
        num_conditions = records_df['Condition'].nunique()
        num_instances = records_df['Instance'].nunique()
        num_timebins = timebin_interval.num_of_bins()
        X_factor_base = X_factor.reshape(num_factors, num_conditions, num_instances, num_timebins)
        records_base = list(
            product(list(range(num_factors)), list(records_df['Condition'].unique()), list(range(num_instances))))
        records_base_df = pd.DataFrame(records_base, columns=['Factor', 'Condition', 'Instance'])
        records_base_df['Timeseries'] = list(X_factor_base.reshape(len(records_base_df), num_timebins))
        return records_base_df
