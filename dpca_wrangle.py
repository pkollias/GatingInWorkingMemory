import sys
from rec_analyses import *


def main():

    args_version = sys.argv[1:]

    """ factor=, fr=, counts_thr=, fr_thr=, [overwrite=] """
    # args_version = ['factor=StimulusGatingPreBool', 'fr=ConcatFactor2', 'counts_thr=15', 'fr_thr=100']
    # args_version = ['job_id=0', 'overwrite=True']

    # load analysis parameters
    version = job_scheduler(args_version, args_from_parse_func)

    # create analysis object
    dpca = DemixedPrincipalComponent(DataBase(['trials', 'units', 'events', 'conditions']), version)
    db, md = dpca.db, dpca.db.md

    # overwrite check
    target_filename = dpca.get_path_base('valid_units', dpca.get_wrangle_stem())
    print(target_filename)
    if path.exists(target_filename) and ('overwrite' not in version.keys() or not eval(version['overwrite'])):
        exit()

    # mark units as valid based on number of events and
    units = db.tables['units']
    valid_units_behavioral_lists = units.apply(lambda row: dpca.assess_unit_events(row.name), axis=1)

    md.np_saver(valid_units_behavioral_lists, target_filename)


def args_from_parse_func(parse_version):

    args_version_list = []

    # for factor, counts_thr in [('GatedStimulusPostDistMemory', '6'), ('GatedStimulusPostDistMemory', '8'),
    #                            ('GatedStimulusPostDistMemory', '10'), ('GatedStimulusPostDistMemory', '12')]:
    #     args_factor = ['factor={0:s}'.format(factor)]
    #     args_fr = ['fr=ConcatFactor2']
    #     args_counts_thr = ['counts_thr={0:s}'.format(counts_thr)]
    #     args_fr_thr = ['fr_thr=100']
    #     args_version_list.extend(list(map(list, list(product(args_factor, args_fr, args_counts_thr, args_fr_thr)))))

    for counts_thr in ['10', '12', '16']:
        # for factor in ['GatPostStimulusRuleStim', 'TargPostStimulusRuleStim',
        #                'GatingPreBool', 'StimulusGating', 'StimulusGatingPreBool']:
        for factor in ['StimulusGatingNoTarget']:
            args_factor = ['factor={0:s}'.format(factor)]
            args_fr = ['fr=ConcatFactor2']
            args_counts_thr = ['counts_thr={0:s}'.format(counts_thr)]
            args_fr_thr = ['fr_thr=100']
            args_version_list.extend(list(map(list, list(product(args_factor, args_fr, args_counts_thr, args_fr_thr)))))

    args_version_from_job = args_version_list[int(parse_version['job_id'])]
    if 'overwrite' in parse_version.keys():
        args_version_from_job.append('overwrite={0:s}'.format(parse_version['overwrite']))

    return args_version_from_job


main()
