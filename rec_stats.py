import copy
from metadata import *

class PopulationBehavioralTimeseries():

    def __init__(self, shape=(0, 0, 0, 0), dim_order=['Units', 'Conditions', 'Instances', 'Timebins'], condition_levels=[], condition_labels=[]):

        self.data = np.empty(shape)
        self.base_shape = self.data.shape
        self.dim_order = dim_order
        self.format = 'Base'

        self.unit_inds = []
        self.condition_levels = condition_levels
        self.condition_labels = condition_labels
        self.timebins = {'version_fr': np.nan,
                         'timebin': np.nan,
                         'timestep': np.nan}

    def get_current_shape(self):

        return self.data.shape

    def get_base_shape(self):

        return self.base_shape

    def get_data(self):

        return self.data

    def set_data(self, data):

        self.data = data
        self.base_shape = self.data.shape

    def set_no_base_data(self, data):
        self.data = data

    def get_dim_order(self):

        return self.dim_order

    def set_dim_order(self, dim_order):

        self.dim_order = dim_order

    def get_unit_inds(self):

        return self.unit_inds

    def set_unit_inds(self, unit_inds):

        self.unit_inds = unit_inds

    def get_condition_levels(self):

        return self.condition_levels

    def set_condition_levels(self, condition_levels):

        self.condition_levels = condition_levels

    def get_condition_labels(self):

        return self.condition_labels

    def set_condition_labels(self, condition_labels):

        self.condition_labels = condition_labels

    def derive_unit(self, unit_ind=np.nan):
        return UnitBehavioralTimeseries(shape=(0,) + tuple(self.get_current_shape()[2:]), unit_ind=unit_ind)

    def add_unit(self, ubt):
        self.data = np.insert(self.data, self.get_current_shape()[0], ubt.data, axis=0)
        self.unit_inds.append(ubt.get_unit_inds())
        self.base_shape = self.data.shape
        if not bool(self.get_condition_levels()):
            self.set_condition_levels(ubt.get_condition_levels())

    def base_to_PCA(self):

        pbt_pca = copy.deepcopy(self)
        base_shape = self.get_base_shape()
        dim_order = self.get_dim_order()
        pbt_pca.set_no_base_data(self.data.reshape(base_shape[0], np.product(base_shape[1:])).transpose())
        pbt_pca.set_dim_order(['_'.join(dim_order[1:]), dim_order[0]])
        pbt_pca.format = 'PCA'

        return pbt_pca

    def PCA_to_base(self):

        pbt = copy.deepcopy(self)
        dim_order = self.get_dim_order()
        pbt.set_no_base_data(self.data.transpose().reshape(self.get_base_shape()))
        pbt.set_dim_order([dim_order[-1]] + dim_order[:-1][0].split('_'))
        pbt.format = 'Base'

        return pbt


class UnitBehavioralTimeseries(PopulationBehavioralTimeseries):

    def __init__(self, shape=(0, 0, 0), unit_ind=np.nan, condition_levels=[]):
        super(UnitBehavioralTimeseries, self).__init__(shape=((1,) + shape), condition_levels=condition_levels)
        self.set_unit_inds(unit_ind)

    def set_data(self, data):

        self.data = data.reshape((1, ) + data.shape)
        self.base_shape = self.data.shape

    def derive_unit_condition(self, condition_level):

        return UnitConditionTimeseries(shape=(0,) + tuple(self.get_current_shape()[3:]), condition_level=condition_level)

    def add_condition(self, uct):

        self.data = np.insert(self.data, self.get_current_shape()[1], uct.data, axis=1)
        self.condition_levels.append(uct.get_condition_levels())
        self.base_shape = self.data.shape





class UnitConditionTimeseries(UnitBehavioralTimeseries):

    def __init__(self, shape=(0, 0), condition_level=np.nan):
        super(UnitConditionTimeseries, self).__init__(shape=(1,) + shape)
        self.set_condition_levels(condition_level)

    def set_data(self, data):

        self.data = data.reshape((1, 1,) + data.shape)
        self.base_shape = self.data.shape
