from rec_analyses import *
import matplotlib.pyplot as plt



def unit_wv_plot(u_iloc, ax=None):

    # u_iloc = 0

    md = MetaData()
    db = md.db_base_loader(['units'])
    units = db['units']

    sess, channum, unitnum = units.iloc[u_iloc].name
    wv = md.np_loader(md.spike_dest_path('wvs', sess, channum, unitnum))
    wv_mean = np.mean(wv, axis=0)

    if ax is None:
        fig, ax = plt.subplots()
    ax.plot(wv_mean)

    plt.box(False)
