import sys
from rec_analyses import *
from sklearn.pipeline import Pipeline
from rec_utils import fbt_df_from_PCA

def main():
    """ class_list, balance_list, fr, counts_thr, area_list, subject, area, mode, mode_seed, pca_mode, [overwrite] """
    args_version = sys.argv[1:]
    # args_version = ['job_id=0', 'overwrite=True']
    version = job_scheduler(args_version, args_from_parse_func)

    class_list = version['class_list'].split('_')
    balance_list = version['balance_list'].split('_')

    def version_modifier(class_i, balance_i):
        version_mod = version.copy()
        version_mod['class'] = class_i
        version_mod['balance'] = balance_i
        return version_mod

    # create analysis object
    classifier_list = [ClassificationAnalysis(DataBase([]), version_modifier(class_i, balance_i))
                       for class_i, balance_i
                       in zip(class_list, balance_list)]
    md = MetaData()

    # overwrite check
    stem = classifier_list[0].get_assemble_stem()
    fname = 'pca' + ('_mean' if version['pca_mode'] == 'mean' else '')
    target_filename = classifier_list[0].get_path_base(fname, stem, cross=True)
    print(target_filename, file=sys.stdout)
    if path.exists(target_filename) and ('overwrite' not in version.keys() or not eval(version['overwrite'])):
        exit()

    # load pbt and crop timeseries
    pbt_list = [md.np_loader(class_i.get_path_base('pbt', class_i.get_assemble_stem())) for class_i in classifier_list]
    for pbt_i in pbt_list:
        pbt_i.crop_timeseries(-50, 1000)

    # conduct pca on combined pbts
    pbt = pbt_list[0].init_with_df(pd.concat([pbt_i.df for pbt_i in pbt_list]))
    X, records = pbt.average_instances(['Unit', 'Unit_Code', 'Condition']).to_PCA_array() if version['pca_mode'] == 'mean' else pbt.to_PCA_array()
    scaler = StandardScaler()
    X_sc = scaler.fit_transform(X)
    pca = PCA()
    pca.fit(X_sc)

    # save pca pipeline
    pipeline = Pipeline([('scaler', scaler), ('pca', pca)])
    md.np_saver(pipeline, target_filename)

    # apply pca to individual pbts and save
    for classifier_i, pbt_i in zip(classifier_list, pbt_list):
        X_i, records_i = pbt_i.to_PCA_array()
        X_factor_i = pipeline.transform(X_i)
        fbt_i = FactorBehavioralTimeseries(fbt_df_from_PCA(X_factor_i, records_i,
                                                           X_factor_i.shape[1], pbt_i.timebin_interval),
                                           pbt_i.condition_labels, pbt_i.timebin_interval)
        stem = tuple(list(classifier_i.get_assemble_stem()) + \
                     ['_'.join([classifier_i.version[key] for key in ['class', 'balance']])])
        fname = 'fbt' + ('_mean' if version['pca_mode'] == 'mean' else '')
        target_filename = classifier_i.get_path_base(fname, stem, cross=True)
        md.np_saver(fbt_i, target_filename)


def args_from_parse_func(parse_version):

    args_version_list = []

    for args_class_list, args_balance_list in [(['class_list=Stimulus_GatedStimulus'], ['balance_list=StageGatingPrePostMemory_StageGatingPrePostSensory'])]:
        for counts_thr in ['12']:  # ['15', '12', '9']:
            for area_list in ['PFC', 'Stri', 'IT']:
                for area in area_list.split('_'):
                    args_fr = ['fr=ConcatFactor2']
                    args_counts_thr = ['counts_thr={0:s}'.format(counts_thr)]
                    args_area_list = ['area_list={0:s}'.format(area_list)]
                    args_subject = ['subject=Gonzo_Oscar']
                    args_area = ['area={0:s}'.format(area)]
                    args_mode = ['mode=Normal']
                    args_mode_seed = ['mode_seed=0']
                    args_pca_mode = ['pca_mode=mean']
                    args_version_list.extend(list(map(list, list(product(args_class_list, args_balance_list, args_fr, args_counts_thr,
                                                                         args_area_list, args_subject, args_area,
                                                                         args_mode, args_mode_seed, args_pca_mode)))))

    args_version_from_job = args_version_list[int(parse_version['job_id'])]
    if 'overwrite' in parse_version.keys():
        args_version_from_job.append('overwrite={0:s}'.format(parse_version['overwrite']))

    return args_version_from_job


main()
