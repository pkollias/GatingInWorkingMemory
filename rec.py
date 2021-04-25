from __future__ import annotations
import copy
import scipy.sparse
import astropy
from astropy import units as qu
from typing import Callable
from metadata import *


class SamplingMethod:

    def __init__(self, srate: astropy.units.quantity.Quantity=30e3*qu.Hz):

        self.srate = srate.to('Hz')

    def stamp_from_t(self, t: astropy.units.quantity.Quantity):

        return int(t.to('s').value * self.srate.to('Hz').value)

    def t_from_stamp(self, stamp) -> astropy.units.quantity.Quantity:

        return (stamp / self.srate.to('Hz')).to('s')

    def reset_sampling(self, srate: astropy.units.quantity.Quantity):

        self.srate = srate.to('Hz')

    def sampling_point(self, t: astropy.units.quantity.Quantity=None, stamp=None) -> SamplingPoint:

        return SamplingPoint(self, t, stamp)

    def zero(self) -> SamplingPoint:

        return SamplingPoint(self, t=0*qu.s)

    def __str__(self):

        return 'SR:{0}'.format(self.srate)

    def __repr__(self):

        return '{0} at {1}\n{2}'.format(type(self), hex(id(self)), str(self))



class SamplingPoint:

    def __init__(self, sampling: SamplingMethod, t: astropy.units.quantity.Quantity=None, stamp=None):

        self._sampling = sampling

        if t is not None and stamp is None:
            self.t = t.to('s')
            self.stamp = self._sampling.stamp_from_t(t)
        elif stamp is not None and t is None:
            self.stamp = int(stamp)
            self.t = self._sampling.t_from_stamp(stamp)
        elif t is None and stamp is None:
            self.t = None
            self.stamp = None
        else:
            self.t = t.to('s')
            self.stamp = int(stamp)

    def copy(self) -> SamplingPoint:

        return type(self)(self._sampling, self.t, self.stamp)

    def get_offset(self, sp: SamplingPoint) -> SamplingPoint:

        return type(self)(self._sampling, self.t+sp.t, self.stamp+sp.stamp)

    def slide_by_offset(self, sp: SamplingPoint) -> SamplingPoint:

        self.t = self.t + sp.t
        self.stamp = self.stamp + sp.stamp

    def __str__(self):

        return '{0} {1} :{2}:]'.format(str(self._sampling), self.t, self.stamp)

    def __repr__(self):

        return '{0} at {1}\n{2}'.format(type(self), hex(id(self)), str(self))



class SamplingInterval:

    def __init__(self, sp_start: SamplingPoint, sp_end: SamplingPoint):

        self.start = sp_start
        self.end = sp_end

    def copy(self) -> SamplingInterval:

        return type(self)(self.start.copy(), self.end.copy())

    def get_stamp_range_list(self):

        return list(range(self.start.stamp, self.end.stamp))

    def get_t_range(self):

        return [self.start._sampling.t_from_stamp(stamp) for stamp in list(range(self.start.stamp, self.end.stamp))]

    def get_offset(self, value) -> SamplingInterval:

        offset = resolve_to_interval(value)
        return type(self)(self.start.get_offset(offset.start), self.end.get_offset(offset.end))

    def resize_by_offset(self, value) -> SamplingInterval:

        offset = resolve_to_interval(value)
        self.start.slide_by_offset(offset.start)
        self.end.slide_by_offset(offset.end)

    def slice_by_index(self, slice_range: SamplingInterval) -> SamplingInterval:

        return resolve_to_interval(self.start.copy()).get_offset(slice_range)

    def index_of_slice(self, slice_interval: SamplingInterval) -> SamplingInterval:

        return slice_interval.get_offset(SamplingPoint(self.start._sampling, stamp=-self.start.stamp))

    def bins_by_step(self, t_bin: SamplingPoint, t_step: SamplingPoint) -> [SamplingInterval]:

        nbins = round(((self.duration().value - t_bin.t.value + t_step.t.value) / t_step.t.value))
        sm = self.start._sampling
        step_stamps = [int(stamp) for stamp in np.linspace(0, (nbins - 1) * t_step.stamp, nbins)]
        return [self.slice_by_index(SamplingInterval(SamplingPoint(sm, stamp=step_stamps[bin_i]),
                                                     SamplingPoint(sm, stamp=step_stamps[bin_i] + t_bin.stamp)))
                for bin_i in range(nbins)]

    def duration(self) -> astropy.units.quantity.Quantity:

        return self.end.t - self.start.t

    def num_samples(self):

        return self.end.stamp - self.start.stamp

    def __str__(self):

        return '[{0} - {1}]'.format(str(self.start), str(self.end))

    def __repr__(self):

        return '{0} at {1}\n{2}'.format(type(self), hex(id(self)), str(self))



class Session:

    def __init__(self, entry, sampling: SamplingMethod=SamplingMethod(MetaData.base_rate * qu.Hz)):

        self.session_entry = entry
        self._session_sampling = sampling
        session_start = SamplingPoint(self._session_sampling, stamp=0)
        session_end = SamplingPoint(self._session_sampling, stamp=entry.TimeStamps - 1)
        self.session_interval = SamplingInterval(session_start, session_end)

    def __str__(self):

        return self.session_entry.Session

    def __repr__(self):

        return str(self)



class Signal:

    padding_value = np.nan

    def __init__(self, session, sampling: SamplingMethod, data, interval: SamplingInterval, units: astropy.units.core.Unit=None, src_data=None):

        self._session = session
        self._sampling = sampling

        self.data = data
        self.units = units
        self.interval = interval
        self.interval_dict = {}
        # src_data has information of original signal, could be in different sampling format or different slicing,
        # or could be a collection of signals (list or other structure)
        self.src_data = src_data

    def copy(self) -> Signal:

        return type(self)(self._session, self._sampling, copy.deepcopy(self.data), self.interval.copy(), self.units, self.src_data)

    def slice_by_index(self, slice_range: SamplingInterval, keep_src=False) -> Signal:

        return Signal(self._session,
                      self._sampling,
                      copy.deepcopy(self.data[slice_range.start.stamp:slice_range.end.stamp]),
                      self.interval.slice_by_index(slice_range),
                      self.units,
                      (self.src_data if keep_src else None))

    def index_of_slice(self, slice_interval: SamplingInterval) -> SamplingInterval:

        return self.interval.index_of_slice(slice_interval)

    def slice_by_sampling_interval(self, slice_interval: SamplingInterval, keep_src=False) -> Signal:

        return self.slice_by_index(self.index_of_slice(slice_interval), keep_src)



class SignalAggregation():

    def __init__(self):

        self.data_list = []
        self.signal_list = []

    def append_signal(self, signal: Signal) -> SignalAggregation:

        self.data_list.append(signal.data)
        self.signal_list.append(signal)



class SpikeTrain(Signal):

    def __init__(self, session, sampling: SamplingMethod, data: scipy.sparse.coo.coo_matrix, interval: SamplingInterval, src_data=None):

        super(SpikeTrain, self).__init__(session, sampling, data, interval, None, src_data)

    def select_by_index(self, select_range: SamplingInterval, keep_src=False) -> SpikeTrain:

        ind_range = np.searchsorted(self.data.col, [select_range.start.stamp, select_range.end.stamp], 'left')
        data_slice_stamps = self.data.col[ind_range[0]:ind_range[1]]
        data_slice_vals = self.data.data[ind_range[0]:ind_range[1]]
        return SpikeTrain(self._session,
                          self._sampling,
                          stamps_to_sparse(data_slice_stamps, data_slice_vals, self.interval.num_samples()),
                          self.interval,
                          (self.src_data if keep_src else None))

    def slice_by_index(self, slice_range: SamplingInterval, keep_src=False) -> SpikeTrain:

        selection = self.select_by_index(slice_range, keep_src)
        selection_stamps = selection.data.col - slice_range.start.stamp
        selection_vals = selection.data.data
        selection_num_data_points = slice_range.num_samples()
        selection.data = stamps_to_sparse(selection_stamps, selection_vals, selection_num_data_points)
        selection.interval = self.interval.slice_by_index(slice_range)
        return selection

    def firing_rate(self) -> float:

        numspikes = len(self.data.col)
        event_duration = self.interval.duration().value
        return numspikes / event_duration

    def to_MultiSpikeTrain(self) -> MultiSpikeTrain:

        multi_spikedata_vals = [1 for _ in np.arange(len(self.data.data))]
        multi_spikedata = stamps_to_sparse(self.data.col, multi_spikedata_vals, self.interval.num_samples())
        return MultiSpikeTrain(self._session, self._sampling, multi_spikedata, self.interval, qu.Hz, self.src_data)



class MultiSpikeTrain(SpikeTrain):

    def __init__(self, session, sampling: SamplingMethod, data, interval: SamplingInterval, units: astropy.units.core.Unit, src_data=None):

        super(MultiSpikeTrain, self).__init__(session, sampling, data, interval, src_data)
        self.units = units



class SignalSmoothing:

    def __init__(self, func: Callable=signal.convolve, window: np.ndarray=MetaData.filt_win_gauss):

        self.func = func
        self.window = window

    def smoothen(self, data: np.ndarray):

         return self.func(data, self.window, mode='same')/sum(self.window)



class TimebinInterval:

    def __init__(self, timebin, timestep, t_start, t_end):

        self.timebin = timebin
        self.timestep = timestep
        self.t_start = t_start
        self.t_end = t_end

    def num_of_bins(self):

        return int(((self.t_end - self.t_start - self.timebin) / self.timestep) + 1)


    def split_to_bins_onset(self):

        timebin, timestep, t_start, t_end = self.timebin, self.timestep, self.t_start, self.t_end
        return [int(onset)
                for onset
                in np.linspace(t_start, t_end - timebin, self.num_of_bins(), endpoint=True)]

    def split_to_bins_offset(self):

        timebin, timestep, t_start, t_end = self.timebin, self.timestep, self.t_start, self.t_end
        return [int(offset)
                for offset
                in np.linspace(t_start + timebin, t_end, self.num_of_bins(), endpoint=True)]

    def sub_interval(self, t_start, t_end):
        return TimebinInterval(self.timebin, self.timestep, t_start, t_end)


# TODO: implement later, bin, stepped timeseries
class Timeseries:

    def __init__(self):

        pass







### Misc ###

def resolve_to_interval(value):

    if type(value) == SamplingPoint:
        return SamplingInterval(value, value)
    elif type(value) == SamplingInterval:
        return value



def stamps_to_sparse(stamps, vals, num_data_points):

    return scipy.sparse.coo_matrix((vals, ([0 for _ in np.arange(len(stamps))], stamps)), shape=(1, num_data_points))