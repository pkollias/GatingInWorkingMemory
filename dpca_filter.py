import sys
from rec_analyses import *


def main():
    """ factor, fr, counts_thr, area_list, subject, [overwrite] """
    args_version = sys.argv[1:]
    # args_version = ['job_id=0', 'overwrite=True']
    version = job_scheduler(args_version, args_from_parse_func)

    # create analysis object
    dpca = DemixedPrincipalComponent(DataBase(['sessions', 'units', 'events', 'conditions']), version)
    db, md = dpca.db, dpca.db.md

    # overwrite check
    target_filename = dpca.get_path_base('filter', dpca.get_filter_stem())
    print(target_filename)
    if path.exists(target_filename) and ('overwrite' not in version.keys() or not eval(version['overwrite'])):
        exit()

    # create unit selection filters
    units = db.tables['units']
    # valid units
    valid_units_behavioral_lists = md.np_loader(dpca.get_path_base('valid_units', dpca.get_wrangle_stem()))
    valid_units = valid_units_behavioral_lists.apply(bool)
    # single units
    single_units = units['UnitNum'].ne(0) & units['RatingCode'].ne(7)
    # area list units
    area_list_units = units['Area'].isin(dpca.version['area_list'].split('_'))
    # subject units
    sessions = db.tables['sessions']
    subject_sessions = sessions.loc[sessions['Subject'].isin(dpca.version['subject'].split('_'))].index
    subject_units = units['Session'].isin(subject_sessions)
    # apply filters
    db.tables['units'] = units.loc[valid_units & single_units & area_list_units & subject_units]

    # create behavioral_units table
    units_behavioral_lists = valid_units_behavioral_lists.loc[db.tables['units'].index]
    behavioral_units = pd.DataFrame([(sess, channum, unitnum, trialnum, stageindex)
                                     for (sess, channum, unitnum), events_list in units_behavioral_lists.iteritems()
                                     for _, trialnum, stageindex in events_list],
                                    columns=md.proc_imports['behavioral_units']['index'])

    md.np_saver(behavioral_units, target_filename)


def args_from_parse_func(parse_version):

    args_version_list = []

    # for area_list, area in [('PFC', 'PFC'), ('Stri', 'Stri'), ('IT', 'IT')]:
    #     for subject in ['Gonzo', 'Oscar', 'Gonzo_Oscar']:
    #         for factor, counts_thr in [('GatedStimulusPostDistMemory', '6'), ('GatedStimulusPostDistMemory', '8'),
    #                                    ('GatedStimulusPostDistMemory', '10'), ('GatedStimulusPostDistMemory', '12')]:
    #             args_factor = ['factor={0:s}'.format(factor)]
    #             args_fr = ['fr=ConcatFactor2']
    #             args_counts_thr = ['counts_thr={0:s}'.format(counts_thr)]
    #             args_area_list = ['area_list={0:s}'.format(area_list)]
    #             args_subject = ['subject={0:s}'.format(subject)]
    #             args_version_list.extend(list(map(list, list(product(args_factor, args_fr, args_counts_thr, args_area_list, args_subject)))))

    for area_list, area in [('PFC', 'PFC'), ('Stri', 'Stri'), ('IT', 'IT')]:
        for subject in ['Gonzo', 'Oscar', 'Gonzo_Oscar']:
            for counts_thr in ['10', '12', '16']:
                # for factor in ['GatPostStimulusRuleStim', 'TargPostStimulusRuleStim',
                #                'GatingPreBool', 'StimulusGating', 'StimulusGatingPreBool']:
                for factor in ['StimulusGatingNoTarget']:
                    args_factor = ['factor={0:s}'.format(factor)]
                    args_fr = ['fr=ConcatFactor2']
                    args_counts_thr = ['counts_thr={0:s}'.format(counts_thr)]
                    args_area_list = ['area_list={0:s}'.format(area_list)]
                    args_subject = ['subject={0:s}'.format(subject)]
                    args_version_list.extend(list(map(list, list(product(args_factor, args_fr, args_counts_thr, args_area_list, args_subject)))))

    args_version_from_job = args_version_list[int(parse_version['job_id'])]
    if 'overwrite' in parse_version.keys():
        args_version_from_job.append('overwrite={0:s}'.format(parse_version['overwrite']))

    return args_version_from_job


main()
