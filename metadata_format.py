from numpy import linspace



def behunit_params_str(version, timebin, timestep, t_start, t_end):

    return '{0:s}_bin{1:04d}_step{2:04d}_t{3:04d}_t{4:04d}'.format(version, timebin, timestep, t_start, t_end).replace('-', 'n')



def interval_split_to_bins_onset(t_start, t_end, timebin, timestep):

    return [int(onset)
            for onset
            in linspace(t_start,
                        t_end - timebin,
                        int(((t_end - t_start - timebin) / timestep) + 1),
                        endpoint=True)]