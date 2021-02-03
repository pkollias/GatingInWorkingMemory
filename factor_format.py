from metadata_format import *



def factor_filter_params_str(filter_params_list_str, counts_thr, area_list_str):

    return '{0:s}_{1:s}_{2:03d}'.format(filter_params_list_str, area_list_str, counts_thr)