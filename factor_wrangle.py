import sys
from rec import *
from rec_format import *


def main():


    # load analysis parameters
    args = sys.argv
    version_factor = args[1]
    version_fr = args[2]


    # version parameters
    v_factor_params = factor_generate_conditions(version_factor)
    condition_columns = v_factor_params['condition_columns']
    condition_list = v_factor_params['condition_list']

    v_fr_params = anova_version_fr_params(version_fr)
    t_start = v_fr_params['t_start']
    t_end = v_fr_params['t_end']
    timebin = v_fr_params['timebin']
    timestep = v_fr_params['timestep']



    # load data and init vars, tables, and slices
    md = MetaData()
    db = md.db_base_loader(['events', 'units', 'conditions'])
    events, units, conditions = db['events'], db['units'], db['conditions']

    target_filename = md.proc_dest_path(path.join('BehavioralUnits', 'Factorization', version_factor,
                                                  behunit_params_str(version_fr, timebin, timestep, t_start, t_end), 'wrangle'),
                                        'conditions_events_dict.pkl')
    print(target_filename)
    if path.exists(target_filename):
        exit()

    events_index = md.preproc_imports['events']['index']
    events_conditions = pd.merge(events.reset_index(drop=True), conditions, on=events_index).set_index(events_index, drop=False)


    factor_condition_events_dict = {}

    # for every unit
    for ind, unit_index in enumerate(units.index):

        print(ind, unit_index)

        sess, channum, unitnum = unit_index

        # load unit's firing rate data
        src_filename = md.proc_dest_path(path.join('BehavioralUnits', 'FiringRates',
                                                   behunit_params_str(version_fr, timebin, timestep, t_start, t_end)),
                                         'behunit_FR_{0:s}_chan{1:03d}_unit{2:03d}.pkl'.format(sess, channum, unitnum))
        timebin_fr_dict = md.np_loader(src_filename)

        if bool(timebin_fr_dict):
            # get list of event indices for that unit's activity
            event_index_list = list(timebin_fr_dict)

            events_conditions_slice = events_conditions.loc[event_index_list]
            # create column for factor condition values
            events_conditions_slice['Factor_Condition'] = events_conditions_slice.apply(lambda row: tuple([row[col] for col in condition_columns]), axis=1)

            # filter only rows that belong in condition list
            events_conditions_collection = events_conditions_slice[events_conditions_slice['Factor_Condition'].isin(condition_list)]
            # group data by factor condition
            condition_groupping = events_conditions_collection.groupby(events_conditions_collection['Factor_Condition'])

            # dictionary of unit_dicts -> dictionary of factor_conditions_dicts -> dictionary of events_lists/timeseries -> list of firing rates
            # aka timeseries for every event for every condition for every unit
            unit_cond_dict = {cond: event_index_vals
                              for cond, event_index_vals
                              in condition_groupping.groups.items()}
            # whether all conditions in analysis have at least one event
            valid_conditions = not bool(set.difference(set(condition_list), set(unit_cond_dict.keys())))

            unit_dict = {'unit_cond_dict': unit_cond_dict,
                         'valid_conditions': valid_conditions}
        else:
            unit_dict = {'unit_cond_dict': {},
                         'valid_conditions': False}

        factor_condition_events_dict[unit_index] = unit_dict

    md.np_saver(factor_condition_events_dict, target_filename)


main()
