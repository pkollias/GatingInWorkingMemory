def behunit_params_str(version, timebin, timestep, t_start, t_end):

    return '{0:s}_bin{1:04d}_step{2:04d}_t{3:04d}_t{4:04d}'.format(version, timebin, timestep, t_start, t_end).replace('-', 'n')


#######
# Anova

def interaction_term(x_a, x_b):
    x_list = [x_a, x_b]
    x_list.sort()
    return ':'.join(x_list)


def y_bin_strnum_to_str(y):

    return 'bin_{0:04d}'.format(int(y)).replace('-', 'n')


def y_bin_str_to_strnum(y):

    return y.replace('bin_', '').replace('n', '-')


def shuffle_to_name(shuffle_i):
    if shuffle_i == 0:
        return 'observed'
    else:
        return 'shuffle_{0:04d}'.format(shuffle_i - 1)


#######
# Factor

def factor_filter_params_str(filter_params_list_str, counts_thr, area_list_str):

    return '{0:s}_{1:s}_{2:03d}'.format(filter_params_list_str, area_list_str, counts_thr)