# fourier.py
# Module to perfor Fourier analysis of simulation data
# Copyright 2015 Giuseppe Venturini
#
# This file is part of the ahkab simulator.
#
# Ahkab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 2 of the License.
#
# Ahkab is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License v2
# along with ahkab.  If not, see <http://www.gnu.org/licenses/>.

"""
This module offers the functions needed to perform a Fourier analysis of the
results of a simulation.

Module reference
''''''''''''''''

"""
from __future__ import (unicode_literals, absolute_import,
                        division, print_function)

import numpy as np
import numpy.fft as fft

from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.signal.windows import (bartlett, hann, hamming, blackman,
                          blackmanharris, gaussian, kaiser)

from . import options
from . import py3compat

# .FFT <output_var> <START = value> <STOP = value> <NP = value>
#+ <FORMAT = keyword> <WINDOW = keyword> <ALFA = value>
#+ <FREQ = value> <FMIN = value> <FMAX = value>
specs = {'fft':{'tokens':({
                          'label':None,
                          'pos':0,
                          'type':str,
                          'needed':True,
                          'dest':'label',
                          'default':None
                         },
                         {
                          'label':'start',
                          'pos':None,
                          'type':float,
                          'needed':False,
                          'dest':'stop',
                          'default':None
                         },
                         {
                          'label':'stop',
                          'pos':None,
                          'type':float,
                          'needed':False,
                          'dest':'start',
                          'default':None
                         },
                         {
                          'label':'from',
                          'pos':None,
                          'type':float,
                          'needed':False,
                          'dest':'from',
                          'default':None
                         },
                         {
                          'label':'to',
                          'pos':None,
                          'type':float,
                          'needed':False,
                          'dest':'to',
                          'default':None
                         },
                         {
                          'label':'np',
                          'pos':None,
                          'type':float,
                          'needed':False,
                          'dest':'np',
                          'default':0
                         },
                         {
                          'label':'format',
                          'pos':None,
                          'type':str,
                          'needed':False,
                          'dest':'uformat',
                          'default':'NORM'
                         },
                         {
                          'label':'window',
                          'pos':None,
                          'type':str,
                          'needed':False,
                          'dest':'window',
                          'default':options.RECT_WINDOW
                         },
                         {
                          'label':'alpha',
                          'pos':None,
                          'type':float,
                          'needed':False,
                          'dest':'alpha',
                          'default':3.
                         },
                         {
                          'label':'freq',
                          'pos':None,
                          'type':float,
                          'needed':False,
                          'dest':'freq',
                          'default':None
                         },
                         {
                          'label':'fmin',
                          'pos':None,
                          'type':float,
                          'needed':False,
                          'dest':'fmin',
                          'default':None
                         },
                         {
                          'label':'fmax',
                          'pos':None,
                          'type':float,
                          'needed':False,
                          'dest':'fmax',
                          'default':None
                         }
                        )
               }
           }

#.FOUR <freq> <x1> [<x2> [<x3 ...]]
def fourier(label, tran_results, fund):
    """Fourier analysis of the time evolution of a variable.

    In particular, the function uses the first 10 multiples of the fundamental
    frequency and a rectangular window.

    A variable amount of time data is used, resampled with a fixed
    time step. The length of the data is decided as follows:

    * The data should be taken from the end of the simulation, so that if there
      is any build-up or stabilization process, the Fourier analysis is not
      affected (or less affected) by it.
    * At least 1 period of the fundamental should be used.
    * Not more than 50% of the total simulation time should be used, if possible.
    * Respecting the above, as much data as possible should be used, as it leads
      to more accurate results.

    **Parameters:**

    label : str or tuple of str
        The identifier of a variable. Eg. ``'Vn1'`` or ``'I(VS)'``. If ``r`` is
        your ``tran_solution`` object, calling ``r.keys()`` will give you all
        the possible variable names for your result set.
        If a tuple of two identifiers is provided, the difference of the two,
        in the form ``label[0]-label[1]``, will be used.
    tran_results : tran_solution instance
        The TRAN results containing the time data for the ``'label'`` variable.
    fund : float
        The fundamental frequency, in Hertz.

    **Returns:**

    f : ndarray of floats
        The frequencies correspoding to the ``F`` array below.
    F : ndarray of complex data
        The result of the Fourier transform, including DC.
    THD : float
        The total harmonic distortion. This value, for a meaningful case, should
        be in the range (0, 1).
    """
    # we want a good number of periods. Not too big, to leave space for any
    # transient at the beginning, not too small, to have points in the FFT
    if type(label) in py3compat.string_types:
        data = tran_results[label]
    elif type(label) == tuple:
        if len(label) == 1 or (len(label) >= 2 and label[1] is None):
            data = tran_results[label[0]]
            label = label[0]
        else:
            data = tran_results[label[0]] - tran_results[label[1]]
            label = '%s - %s' % (label[0], label[1])
    ps = []
    for i in range(1, 10, 2):
        ps += [(tran_results.tstop - tran_results.tstart) // (i/fund)]
    ps = [i for i in ps if i > 0]
    nperiods = ps[int(np.ceil(len(ps)-1))]
    # setup start & stop accordingly
    start = tran_results.tstop - nperiods/fund
    stop = tran_results.tstop
    # 10 points, remember?
    sampling = 1./(2*10*fund)
    t = np.linspace(start, stop, (stop-start)//sampling, endpoint=False)
    idata = InterpolatedUnivariateSpline(tran_results.get_x(), data, k=2)
    f = fft.fftfreq(len(t), sampling)[::nperiods]
    F = fft.rfft(idata(t))[:-1:nperiods]
    f = f[:len(f)/2]
    THD = np.sqrt(sum(abs(F[2:])**2))/abs(F[1])
    return f, F, THD

# .FFT <output_var> <START = value> <STOP = value> <NP = value>
#+ <FORMAT = keyword> <WINDOW = keyword> <ALFA = value>
#+ <FREQ = value> <FMIN = value> <FMAX = value>
def spicefft(label, tran_results, freq=None, **args):
    """FFT analysis of the time evolution of a variable.

    This function is a much more flexible and complete version of the
    :func:`ahkab.fourier.fourier` function.

    The function uses a variable amount of time data, resampled with a fixed
    time step.
    The time interval is specified through the ``start`` and ``stop``
    parameters, if they are not set, all the available data is used.

    The function behaves differently whether the parameter ``freq`` is specified
    or not:

    * If the fundamental frequency ``freq`` (:math:`f` in the following) is
      specified, the function will perform an harmonic analysis, considering
      only the DC component and the harmonics of :math:`f` up to the 9th (ie
      :math:`f`, :math:`2f`, :math:`3f` :math:`\\dots` :math:`9f`).
    * If ``freq`` is left unspecified, a standard FFT analysis is performed,
      starting from :math:`f = 0`, to a frequency :math:`f_{max} =
      1/(2T_{TOT}n_p)`, where :math:`T_{TOT}` is the total length of the
      considered data in seconds and :math:`n_p` is the number of points in the
      FTT, set through the ``np`` parameter to this function.

    **Parameters:**

    label : str, or tuple of str
        The identifier of a variable. Eg. ``'Vn1'`` or ``'I(VS)'``. If ``r`` is
        your ``tran_solution`` object, calling ``r.keys()`` will give you all
        the possible variable names for your result set.
        If a tuple of two identifiers is provided, the difference of the two,
        in the form ``label[0]-label[1]``, will be used.
    tran_results : tran_solution instance
        The TRAN results containing the time data for the ``'label'`` variable.
    freq : float, optional
        The fundamental frequency, in Hertz. If it is specified, the output will
        be limited to the harmonics of this frequency. The THD evaluation will
        also be enabled.
    start : float, optional
        The first time instant to be considered for the transient analysis. If
        unspecified, it will be the beginning of the transient simulation.
    from : float, optional
        Alternative specification of the ``start`` parameter.
    stop : float, optional
        Last time instant to be considered for the FFT analysis. If unspecified,
        it will be the end time of the transient simulation.
    to : float, optional
        Alternative specification of the ``stop`` parameter.
    np : integer
        A power of two that specifies how many points should be used when
        computing the FFT. If it is set to a value that is not a power of 2, it
        will be rounded up to the nearest power of 2. It defaults to 1024. 
    window : str, optional
        The windowing type. The following values are available:

        * 'RECT' for a rectangular window, equivalent to no window at all.
        * 'BART', for a Bartlett window.
        * 'HANN', for a Hanning window.
        * 'HAMM' for a Hamming window.
        * 'BLACK' for a Blackman window.
        * 'HARRIS' for a Blackman-Harris window.
        * 'GAUSS' for a Gaussian window.
        * 'KAISER' for a Kaiser-Bessel window.

        The default is the rectangular window.
    alpha : float, optional
        The :math:`\\sigma` for a gaussian window or the :math:`beta` for a
        Kaiser window. Defaults to 3 and is ignored if a window different from
        Gaussian or Kaiser is selected.
    fmin : float, optional
        Suppress all data below this frequency, expressed in Hz. The suppressed
        data is neither returned nor used to compute the THD (if it is computed
        at all). The DC component is always preserved. Defaults to: return and
        use all data.
    fmax : float, optional
        The dual to ``fmin``, discard data above ``fmax`` and also do not use it
        if computing the THD. Defaults to infinity.

    **Returns:**

    f : ndarray of floats
        The frequencies, including the DC.
    F : ndarray of complex data
        The result of the Fourier transform, including DC.
    THD : float
        The total harmonic distortion, if ``freq`` was specified, ``None``
        otherwise.

    .. plot::

        import pylab
        import numpy as np
        import ahkab
        cir = ahkab.Circuit('Test FOUR and FFT')
        mys = ahkab.time_functions.sin(vo=0, va=1, freq=10e3)
        cir.add_resistor('R1', 'n1', cir.gnd, 1e3)
        cir.add_vsource('V1', 'n1', cir.gnd, 1, function=mys)
        tr = ahkab.new_tran(0, 1e-3, 1e-5, x0=None)
        r = ahkab.run(cir, tr)['tran']

        pylab.figure()
        pylab.title('Comparison of different windowing')
        for w in (ahkab.options.RECT_WINDOW, ahkab.options.BART_WINDOW,
                  ahkab.options.HANN_WINDOW, ahkab.options.HAMM_WINDOW,
                  ahkab.options.BLACK_WINDOW, ahkab.options.HARRIS_WINDOW,
                  ahkab.options.KAISER_WINDOW):
            fsp1, Fsp1, THDsp = ahkab.fourier.spicefft('vn1', r, window=w)
            pylab.plot(fsp1, 20*np.log10(abs(Fsp1)),
                       label=ahkab.options.WINDOWS_NAMES[w])
        pylab.xlabel('Frequency [Hz]'); pylab.ylabel('Magnitude [dB]')
        pylab.xlim((0, 100000))
        pylab.legend()

    """
    if type(label) in py3compat.string_types:
        data = tran_results[label]
    elif type(label) == tuple:
        if len(label) == 1 or (len(label) >= 2 and label[1] is None):
            data = tran_results[label[0]]
            label = label[0]
        else:
            data = tran_results[label[0]] - tran_results[label[1]]
            label = '%s - %s' % (label[0], label[1])
    if 'start' not in args and 'from' not in args:
        start = tran_results.tstart
    elif 'start' in args and args['start'] is not None:
        start = args['start']
    elif args['from'] is not None:
        start = args['from']
    else:
        start = tran_results.tstart

    if 'stop' not in args and 'to' not in args:
        stop = tran_results.tstop
    elif 'stop' in args and args['stop'] is not None:
        stop = args['stop']
    elif args['to'] is not None:
        stop = args['to']
    else:
        stop = tran_results.tstop

    if 'np' not in args:
        np2 = 1024
    else:
        # round to the nearest power of two
        np2 = 2 if args['np'] < 2 else 2**int(np.log2(args['np'] - 1) + 1)

    if 'window' in args:
        window_type = args['window'].upper()
        if window_type not in (options.RECT_WINDOW, options.BART_WINDOW,
                               options.HANN_WINDOW, options.HAMM_WINDOW,
                               options.BLACK_WINDOW, options.HARRIS_WINDOW,
                               options.GAUSS_WINDOW, options.KAISER_WINDOW):
            raise ValueError(('fft(): window may be %s, %s, %s, %s, %s, %s, %s'+
                              ' or %s, got %s') % (options.RECT_WINDOW,
                                                   options.BART_WINDOW,
                                                   options.HANN_WINDOW,
                                                   options.HAMM_WINDOW,
                                                   options.BLACK_WINDOW,
                                                   options.HARRIS_WINDOW,
                                                   options.GAUSS_WINDOW,
                                                   options.KAISER_WINDOW,
                                                   window_type))
    else:
        window_type = options.RECT_WINDOW
    alpha = args['alpha'] if 'alpha' in args else 3.
    fmin = args['fmin'] if 'fmin' in args else None
    fmax = args['fmax'] if 'fmax' in args else None

    if freq:
        nperiods = int((stop-start)*freq)
        # adjust start so we have an integer number
        # of periods
        start = stop - nperiods/freq
        sampling = 1./(2*freq*nperiods)
        t = np.linspace(start, stop, ((stop-start)//sampling + 1),
                        endpoint=False)
    else:
        sampling = (stop-start)/np2
        t = np.linspace(start, stop, np2, endpoint=False)
    idata = InterpolatedUnivariateSpline(tran_results.get_x(),
                                         data, k=2)
    window = {options.RECT_WINDOW: lambda x: 1.,
              options.BART_WINDOW: bartlett,
              options.HANN_WINDOW: hann,
              options.HAMM_WINDOW: hamming,
              options.BLACK_WINDOW: blackman,
              options.HARRIS_WINDOW: blackmanharris,
              options.GAUSS_WINDOW: lambda x: gaussian(x, std=alpha),
              options.KAISER_WINDOW: lambda x: kaiser(x, beta=alpha)}
    f = fft.fftfreq(len(t), sampling)
    f = f[:len(f)/2]
    F = fft.rfft(idata(t)*window[window_type](len(t)))[:-1]
    if freq:
        # downsample
        f = f[::nperiods]
        F = F[::nperiods]
    if fmin is not None and fmin > 0:
        if fmin > f.max():
            raise ValueError("fmin=%g Hz but max(f)=%g Hz" % (fmin, f.max()))
        # keep the DC
        F = np.concatenate((np.atleast_1d(F[0]), F[f >= fmin]))
        f = np.concatenate((np.atleast_1d(f[0]), f[f >= fmin]))
    if fmax is not None:
        if fmax < f[1:].min():
            raise ValueError("fmax=%g Hz but min(f)=%g Hz" %
                             (fmax, f[1:].min()))
        F = F[f <= fmax]
        f = f[f <= fmax]
    if freq:
        THD = np.sqrt(sum(abs(F[2:])**2))/abs(F[1])
    else:
        THD = None
    return f, F, THD

