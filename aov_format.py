from metadata_format import *

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


def filter_df_wrapper(df, column, wrapper, arg):

    mask = wrapper(df[column], arg)
    return {'mask': mask,
            'df': df.loc[mask]}
