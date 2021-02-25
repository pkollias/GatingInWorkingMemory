import copy
from metadata import *
from itertools import product
from operator import itemgetter

class PopulationBehavioralTimeseries():

    base_dim_order = ['Units', 'Conditions', 'Instances', 'Timebins']

    def __init__(self, shape=(0, 0, 0, 0), dim_order=base_dim_order,
                 condition_levels=[], condition_labels=[], timebins={'version_fr': np.nan, 'timebin': np.nan, 'timestep': np.nan}):
        self.data = np.empty(shape)
        self.base_shape = self.data.shape
        self.dim_order = dim_order
        self.format = 'Base'
        self.unit_inds = []
        self.condition_levels = condition_levels
        self.condition_labels = condition_labels
        self.timebins = timebins

    def pbt_from_data(self, data, condition_levels=[], condition_labels=[], timebins={'version_fr': np.nan, 'timebin': np.nan, 'timestep': np.nan}, unit_inds=[]):
        ''' creates pbt given specific data array and assigns base_shape based on data array
            result has Base format and data array has to obey base_dim_order order and size'''
        pbt = PopulationBehavioralTimeseries(condition_levels=condition_levels, condition_labels=condition_labels, timebins=timebins)
        pbt.set_base_data(data, unit_inds)
        return pbt

    def get_current_shape(self):
        return self.data.shape

    def set_base_data(self, data, unit_inds=[]):
        ''' sets base data for Base format '''
        self.set_data(data)
        self.base_shape = self.data.shape
        self.unit_inds = unit_inds

    def set_data(self, data):
        ''' sets data without changing base shape of pbt '''
        self.data = data

    def derive_unit(self, unit_ind=np.nan):
        return UnitBehavioralTimeseries(shape=(0,) + tuple(self.get_current_shape()[2:]), unit_ind=unit_ind, condition_levels=[])

    def add_unit(self, ubt):
        self.data = np.insert(self.data, self.get_current_shape()[0], ubt.data, axis=0)
        self.unit_inds.append(ubt.unit_inds)
        self.base_shape = self.data.shape
        if not bool(self.condition_levels):
            self.condition_levels = ubt.condition_levels

    ### PCA ###

    def base_to_PCA(self):
        ''' converts Base format pbt to PCA format (data, dim_order, format) '''
        pbt_pca = copy.deepcopy(self)
        base_shape = self.base_shape
        dim_order = self.dim_order
        pbt_pca.data = self.data.reshape(base_shape[0], np.product(base_shape[1:])).transpose()
        pbt_pca.dim_order = ['_'.join(dim_order[1:]), dim_order[0]]
        pbt_pca.format = 'PCA'
        return pbt_pca

    ###

    def average_instances(self):
        ''' finds Instances dimension and averages across instances '''
        pbt = copy.deepcopy(self)
        base_shape = pbt.base_shape
        new_base_shape = (base_shape[0], base_shape[1], 1, base_shape[3])
        instances_dim_index = pbt.dim_order.index('Instances')
        cur_shape = pbt.get_current_shape()
        new_shape_list = list(cur_shape)
        new_shape_list[instances_dim_index] = 1
        new_shape = tuple(new_shape_list)
        pbt.data = pbt.data.mean(axis=instances_dim_index).reshape(new_shape)
        pbt.base_shape = new_base_shape
        return pbt

    def get_data_without_instance_dim(self):
        ''' returns data without instance dim, averages if multiple instances '''
        pbt_m_inst = self.average_instances()
        instances_dim_index = pbt_m_inst.dim_order.index('Instances')
        new_shape_list = list(pbt_m_inst.get_current_shape())
        new_shape_list.pop(instances_dim_index)
        new_shape = tuple(new_shape_list)
        return pbt_m_inst.data.reshape(new_shape)


    ### Conditions ###

    def base_to_conditions_unfold(self):
        ''' unfold conditions from 1-dimensional tuples into n-dimensional
            allowing array expansion and reshaping for same number of conditions'''
        pbt_cond_nd = copy.deepcopy(self)
        # params and structs for calculations
        num_condition_dims = len(self.condition_labels)
        # turn list of condition levels into list of condition-level tuples (each tuple has all conditions of one dim)
        unzipped_conditions = list(zip(*self.condition_levels))
        # turns unzipped list of tuples into list of sets (drops duplicates)
        condition_levels_by_dim = [list(set(dim_i_levels)) for dim_i_levels in unzipped_conditions]
        # number of level for each set
        num_condition_levels_by_dim = [len(levels) for levels in condition_levels_by_dim]

        # data
        temp_data = self.data.transpose((1, 0, 2, 3))
        if num_condition_dims > 1:
            # initialize new data table with increased number of condition dimensions
            pbt_cond_nd.data = np.empty(tuple(num_condition_levels_by_dim) + temp_data.shape[1:])
            # for every compound condition
            for condition in product(*condition_levels_by_dim):
                # find index in old list
                condition_1d_index = self.condition_levels.index(condition)
                condition_nd_index = [levels.index(condition[dim_index]) for dim_index, levels in
                                      enumerate(condition_levels_by_dim)]
                pbt_cond_nd.data[tuple(condition_nd_index + [slice(None) for _ in range(3)])] = temp_data[
                    condition_1d_index, ...]
        else:
            pbt_cond_nd.data = temp_data

        units_dim_index = num_condition_dims
        instances_dim_index = units_dim_index + 1
        timebin_dim_index = units_dim_index + 2
        pbt_cond_nd.data = pbt_cond_nd.data.transpose(tuple([units_dim_index, instances_dim_index] + list(range(num_condition_dims)) + [timebin_dim_index]))
        # dim_order
        pbt_cond_nd.dim_order = ['Units', 'Instances'] + ['Conditions_' + condition for condition in self.condition_labels] + ['Timebins']
        # condition levels
        pbt_cond_nd.condition_levels = condition_levels_by_dim
        # base
        pbt_cond_nd.format = 'Conditions'
        return pbt_cond_nd

    ### dPCA ###

    def conditions_to_dPCA(self):
        pbt_dpca = copy.deepcopy(self)
        dim_order = pbt_dpca.dim_order
        new_dim_order_index = tuple([1, 0] + list(range(2, len(dim_order))))
        pbt_dpca.data = pbt_dpca.data.transpose(new_dim_order_index)
        pbt_dpca.dim_order = itemgetter(*new_dim_order_index)(pbt_dpca.dim_order)
        pbt_dpca.format = 'dPCA'
        return pbt_dpca

    def base_to_dPCA(self):
        return self.base_to_conditions_unfold().conditions_to_dPCA()


class UnitBehavioralTimeseries(PopulationBehavioralTimeseries):

    def __init__(self, shape=(0, 0, 0), unit_ind=np.nan, condition_levels=[]):
        super(UnitBehavioralTimeseries, self).__init__(shape=((1,) + shape), condition_levels=condition_levels)
        self.unit_inds = unit_ind

    def set_base_data(self, data):
        self.data = data.reshape((1, ) + data.shape)
        self.base_shape = self.data.shape

    def derive_unit_condition(self, condition_levels):
        return UnitConditionTimeseries(shape=(0,) + tuple(self.get_current_shape()[3:]), condition_levels=condition_levels)

    def add_condition(self, uct):
        self.data = np.insert(self.data, self.get_current_shape()[1], uct.data, axis=1)
        self.condition_levels.append(uct.condition_levels)
        self.base_shape = self.data.shape


class UnitConditionTimeseries(UnitBehavioralTimeseries):

    def __init__(self, shape=(0, 0), condition_levels=np.nan):
        super(UnitConditionTimeseries, self).__init__(shape=(1,) + shape)
        self.condition_levels = condition_levels

    def set_base_data(self, data):

        self.data = data.reshape((1, 1,) + data.shape)
        self.base_shape = self.data.shape




class FactorizedBehavioralTimeseries():

    base_dim_order = ['Conditions_Instances', 'Units', 'Timebins']

    def __init__(self, X_factor, condition_levels=[], condition_labels=[], timebins={}):
        self.X_factor = X_factor
        self.condition_levels=condition_levels
        self.condition_labels = condition_labels
        self.timebins = timebins

    def from_PCA(X_pca, base_shape, condition_levels=[], condition_labels=[], timebins={'version_fr': np.nan, 'timebin': np.nan, 'timestep': np.nan}):
        # convert matrix to fbt base dim_order ['Conditions_Instances', 'Factors', 'Timebins']
        # by first converting to pbt (['Units', 'Conditions', 'Instances', 'Timebins']) and then reshaping and transposing
        X_factor = X_pca.transpose().reshape(base_shape).transpose((1, 2, 0, 3)).reshape((base_shape[1] * base_shape[2], base_shape[0], base_shape[3]))
        fbt_condition_levels = [cond for cond in condition_levels for _ in range(base_shape[2])]
        fbt = FactorizedBehavioralTimeseries(X_factor, fbt_condition_levels, condition_labels, timebins)
        return fbt

    def from_dPCA(X_dpca, condition_levels=[], condition_labels=[], timebins={'version_fr': np.nan, 'timebin': np.nan, 'timestep': np.nan}):
        # convert matrix to fbt base dim_order ['Conditions', 'Factors', 'Timebins']
        dpca_shape = X_dpca.shape
        new_shape = tuple([dpca_shape[0], np.product(dpca_shape[1:-1]), dpca_shape[-1]])
        X_factor = X_dpca.reshape(new_shape).transpose((1, 0, 2))
        fbt_condition_levels = list(product(*condition_levels))
        fbt = FactorizedBehavioralTimeseries(X_factor, fbt_condition_levels, condition_labels, timebins)
        return fbt
