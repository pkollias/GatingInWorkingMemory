import sys
from scipy.io import loadmat
import scipy.sparse
from metadata import *

'''
Prereqs
updated sessions and units pkl files

Argument list:
[fsroot, session_id]
'''


def main():

    args = sys.argv
    sess_iloc = int(args[2])

    # INIT
    md = MetaData()
    db = md.db_base_loader(['sessions', 'units'])
    sessions, units = db['sessions'], db['units']

    # import ts and wvs
    session = sessions.iloc[sess_iloc].Session
    subject = sessions.loc[session].Subject
    date = sessions.loc[session].Date
    sessid = sessions.loc[session].SessID
    units_slice = units.loc[session]
    chan_grouper = units_slice.reset_index(drop=True).groupby('ChanNum')
    for channum, chan_group in chan_grouper:

        sortmat = loadmat(md.preproc_dat_spike_path(subject, str(date),
                                                    'gatingwm_{0}_{1}_{2:03d}_chan{3:03d}.mat'.
                                                    format(subject.lower(), date, sessid, channum)))

        for unitnum in chan_group['UnitNum']:

            unit_inds = (sortmat['units'] == unitnum).reshape((-1,))
            ts = np.array(flatten(sortmat['ts'][unit_inds]))
            wvs = sortmat['wvs'][unit_inds, ]
            numdatapoints = sessions.loc[session].TimeStamps
            ts_sparse = scipy.sparse.coo_matrix(([True for _ in range(len(ts))], ([0 for _ in range(len(ts))], (ts-1).T)),
                                                shape=(1, numdatapoints))
            md.np_saver(ts_sparse, md.spike_dest_path('ts', session, channum, unitnum))
            md.np_saver(wvs, md.spike_dest_path('wvs', session, channum, unitnum))


main()
