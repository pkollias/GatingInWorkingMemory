from metadata import *

class PopulationBehavioralTimeseries:

    def __init__(self, shape=(0, 0, 0, 0), dim_order=['Units', 'Conditions', 'Instances', 'Timebins']):

        self.data = np.empty(shape)
        self.dim_order = dim_order

        self.unit_ind_list = []
        self.conditions_list = []
        self.conditions_columns = []
        self.num_instances = np.nan
        self.timebins = {'version_fr': np.nan,
                         'timebin': np.nan,
                         'timestep': np.nan}

