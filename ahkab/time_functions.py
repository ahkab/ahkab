# -*- coding: iso-8859-1 -*-
# time_functions.py
# Time functions for independent sources
# Copyright 2006 Giuseppe Venturini

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
This module contains several basic time functions.

The classes that are found in module are useful to provide a time-varying
characteristic to independent sources.

Notice that the time functions are not restricted to those provided here, the
user is welcome to provide his own.
Implementing a custom time function is easy and common practice, as long as you
are interfacing to the simulator through Python. Please see the dedicated section
:ref:`define-custom-time-functions` below.

Classes defined in this module
------------------------------

.. autosummary::
    pulse
    pwl
    sin
    exp
    sffm
    am

Supplying a time function to an independent source
--------------------------------------------------

Providing a time-dependent characteristic to an independent source is very
simple and probably best explained with an example.

Let's say we wish to define a sinusoidal voltage source with no offset,
amplitude 5V and 1kHz frequency.

It is done in two steps:

* first we define the time function with the built-in class
  :class:`ahkab.time_functions.sin`:

    .. code-block:: python

        sin1k = time_functions.sin(vo=0, va=5, freq=1e3)

* Then we define the voltage source and we assign the time function to it:

    .. code-block:: python

        cir.add_vsource('V1', 'n1', cir.gnd, 1, function=mys)

In the example above, the sine wave is assigned to a voltage source ``'V1'``,
that gets added to a circuit ``cir`` (not shown).

.. _define-custom-time-functions:

Defining custom time functions
------------------------------

Defining a custom time function is easy, all you need is either:

* A function that takes a ``float`` (the time) and returns the function
  value,
* An instance with a ``__call__(self, time)`` method. This solution
  allows having internal parameters, typically set through the constructor.

In both cases, in time-based simulations, the simulator will call the object at
every time step, supplying a single parameter, the simulation time (``time`` in
the following, of type ``float``).

In turn, the simulator expects to receive as return value a ``float``,
corresponding to the value of the time-dependent function at the time specified
by the ``time`` variable.

If the time-dependent function is used to define the characteristics of a
voltage source (:class:`VSource`), its return value has to be expressed in Volt.
In the case of a current source (:class:`ISource`), the return value is to be
expressed in Ampere.

The standard notation applies.

As an example, we'll define a custom time-dependent voltage source, having a
:math:`\\mathrm{sinc}(ft)` characteristic. In this example, :math:`f` has a
value of 10kHz.

First we define the time function, in this case we'll do that through the Python
``lambda`` construct.

.. code-block:: python

    mys = lambda t: 1 if not t else math.sin(math.pi*1e4*t)/(math.pi*1e4*t)

Then, we define the circuit -- a very simple one in this case -- and assign our
``mys`` function to ``V1``. In the following circuit, we simply apply the
voltage from ``V1`` to a resistor ``R1``.

.. code-block:: python

    import ahkab
    cir = ahkab.Circuit('Test custom time functions')
    cir.add_resistor('R1', 'n1', cir.gnd, 1e3)
    cir.add_vsource('V1', 'n1', cir.gnd, 1, function=mys)
    tr = ahkab.new_tran(0, 1e-3, 1e-5, x0=None)
    r = ahkab.run(cir, tr)['tran']

Plotting ``Vn1`` and the expected result (:math:`\\mathrm{sinc}(ft)`) we
get:

.. plot::

    import math
    import numpy as np
    import pylab
    import ahkab
    cir = ahkab.Circuit('Test custom time functions')
    cir.add_resistor('R1', 'n1', cir.gnd, 1e3)
    mys = lambda t: 1 if not t else math.sin(math.pi*1e4*t)/(math.pi*1e4*t)
    cir.add_vsource('V1', 'n1', cir.gnd, 1, function=mys)
    tr = ahkab.new_tran(0, 1e-3, 1e-5, x0=None)
    r = ahkab.run(cir, tr)['tran']
    t = r.get_x()
    pylab.hold(True)
    pylab.plot(t, r['vn1'], 'o', ms=3, label='V(n1) (simulation)')
    npsin1k = np.frompyfunc(mys, 1, 1)
    pylab.plot(t, npsin1k(t), label='sinc(ft) (theory)')
    pylab.legend()
    pylab.xlabel('t [s]')
    pylab.ylabel('Voltage [V]')


Module reference
----------------
"""

from __future__ import (unicode_literals, absolute_import,
                        division, print_function)

import math

from scipy.interpolate import InterpolatedUnivariateSpline

time_fun_specs = {'sin': { #VO VA FREQ TD THETA
    'tokens': ({
               'label': 'vo',
               'pos': 0,
               'type': float,
               'needed': True,
               'dest': 'vo',
               'default': None
               },
               {
               'label': 'va',
               'pos': 1,
               'type': float,
               'needed': True,
               'dest': 'va',
               'default': None
               },
               {
               'label': 'freq',
               'pos': 2,
               'type': float,
               'needed': True,
               'dest': 'freq',
               'default': None
               },
               {
               'label': 'td',
               'pos': 3,
               'type': float,
               'needed': False,
               'dest': 'td',
               'default': 0.
               },
               {
               'label': 'theta',
               'pos': 4,
               'type': float,
               'needed': False,
               'dest': 'theta',
               'default': 0
               }
               )
        },'exp': { #EXP(V1 V2 TD1 TAU1 TD2 TAU2)
    'tokens': ({
               'label': 'v1',
               'pos': 0,
               'type': float,
               'needed': True,
               'dest': 'v1',
               'default': None
               },
               {
               'label': 'v2',
               'pos': 1,
               'type': float,
               'needed': True,
               'dest': 'v2',
               'default': None
               },
               {
               'label': 'td1',
               'pos': 2,
               'type': float,
               'needed': False,
               'dest': 'td1',
               'default': 0.
               },
               {
               'label': 'tau1',
               'pos': 3,
               'type': float,
               'needed': True,
               'dest': 'tau1',
               'default': None
               },
               {
               'label': 'td2',
               'pos': 4,
               'type': float,
               'needed': False,
               'dest': 'td2',
               'default': float('inf')
               },
               {
               'label': 'tau2',
               'pos': 5,
               'type': float,
               'needed': False,
               'dest': 'tau2',
               'default': float('inf')
               }
               )
        },'pulse': { #PULSE(V1 V2 TD TR TF PW PER)
    'tokens': ({
               'label': 'v1',
               'pos': 0,
               'type': float,
               'needed': True,
               'dest': 'v1',
               'default': None
               },
               {
               'label': 'v2',
               'pos': 1,
               'type': float,
               'needed': True,
               'dest': 'v2',
               'default': None
               },
               {
               'label': 'td',
               'pos': 2,
               'type': float,
               'needed': False,
               'dest': 'td',
               'default': 0.
               },
               {
               'label': 'tr',
               'pos': 3,
               'type': float,
               'needed': True,
               'dest': 'tr',
               'default': None
               },
               {
               'label': 'tf',
               'pos': 4,
               'type': float,
               'needed': True,
               'dest': 'tf',
               'default': None
               },
               {
               'label': 'pw',
               'pos': 5,
               'type': float,
               'needed': True,
               'dest': 'pw',
               'default': None
               },
               {
               'label': 'per',
               'pos': 6,
               'type': float,
               'needed': True,
               'dest': 'per',
               'default': None
               })
        }, 'sffm': { ## SFFM(VO VA FC MDI FS TD)
    'tokens': ({
               'label': 'vo',
               'pos': 0,
               'type': float,
               'needed': True,
               'dest': 'vo',
               'default': None
               },
               {
               'label': 'va',
               'pos': 1,
               'type': float,
               'needed': True,
               'dest': 'va',
               'default': None
               },
               {
               'label': 'fc',
               'pos': 2,
               'type': float,
               'needed': False,
               'dest': 'fc',
               'default': None
               },
               {
               'label': 'mdi',
               'pos': 3,
               'type': float,
               'needed': True,
               'dest': 'mdi',
               'default': None
               },
               {
               'label': 'fs',
               'pos': 4,
               'type': float,
               'needed': True,
               'dest': 'fs',
               'default': None
               },
               {
               'label': 'td',
               'pos': 5,
               'type': float,
               'needed': False,
               'dest': 'td',
               'default': 0.
               })
        }, 'am': { #AM(sa oc fm fc [td])
    'tokens': ({
               'label': 'sa',
               'pos': 0,
               'type': float,
               'needed': True,
               'dest': 'sa',
               'default': None
               },
               {
               'label': 'oc',
               'pos': 1,
               'type': float,
               'needed': True,
               'dest': 'oc',
               'default': None
               },
               {
               'label': 'fm',
               'pos': 2,
               'type': float,
               'needed': True,
               'dest': 'fm',
               'default': None
               },
               {
               'label': 'fc',
               'pos': 3,
               'type': float,
               'needed': True,
               'dest': 'fc',
               'default': None
               },
               {
               'label': 'td',
               'pos': 4,
               'type': float,
               'needed': False,
               'dest': 'td',
               'default': None
               })
        }
}

#
# Functions for time dependent sources  #
#

class pulse(object):
    """Square wave aka pulse function

    .. image:: images/elem/pulse.svg

    **Parameters:**

    v1 : float
        Square wave low value.
    v2 : float
        Square wave high value.
    td : float
        Delay time to the first ramp, in seconds. Negative values are considered
        as zero.
    tr : float
        Rise time in seconds, from the low value ``v1`` to the pulse high value
        ``v2``.
    tf : float
        Fall time in seconds, from the pulse high value ``v2`` to the low value
        ``v1``.
    pw : float
        Pulse width in seconds.
    per : float
        Periodicity interval in seconds.
    """
    # PULSE(V1 V2 TD TR TF PW PER)

    def __init__(self, v1, v2, td, tr, pw, tf, per):
        self.v1 = v1
        self.v2 = v2
        self.td = max(td, 0.0)
        self.per = per
        self.tr = tr
        self.tf = tf
        self.pw = pw
        self._type = "V"

    def __call__(self, time):
        """Evaluate the pulse function at the given time."""
        if time is None:
            time = 0
        time = time - self.per * int(time / self.per)
        if time < self.td:
            return self.v1
        elif time < self.td + self.tr:
            return self.v1 + ((self.v2 - self.v1) / (self.tr)) * (time - self.td)
        elif time < self.td + self.tr + self.pw:
            return self.v2
        elif time < self.td + self.tr + self.pw + self.tf:
            return self.v2 + ((self.v1 - self.v2) / (self.tf)) * (time - (self.td + self.tr + self.pw))
        else:
            return self.v1

    def __str__(self):
        return "type=pulse " + \
            self._type.lower() + "1=" + str(self.v1) + " " + \
            self._type.lower() + "2=" + str(self.v2) + \
            " td=" + str(self.td) + " per=" + str(self.per) + \
            " tr=" + str(self.tr) + " tf=" + str(self.tf) + \
            " pw=" + str(self.pw)


class sin(object):
    """Sine wave

    .. image:: images/elem/sin.svg

    Mathematically, the sine wave function is defined as:

    * :math:`t < t_d`:

    .. math::

        f(t) = v_o + v_a \\sin\\left(\\pi \\phi/180 \\right)

    * :math:`t \\ge t_d`:

    .. math::

        f(t) = v_o + v_a \\exp\\left[-(t - t_d)\,\\theta \\right] \\sin\\left[2 \\pi f (t - t_d) + \\pi \\phi/180\\right]

    **Parameters:**

    vo : float
        Offset value.
    va : float
        Amplitude.
    freq : float
        Sine frequency in Hz.
    td : float, optional
        time delay before beginning the sinusoidal time variation, in seconds. Defaults to 0.
    theta : float optional
        damping factor in 1/s. Defaults to 0 (no damping).
    phi : float, optional
        Phase delay in degrees. Defaults to 0 (no phase delay).

    .. note::

        This implementation is consistent with the SPICE simulator, other simulators use
        different formulae.

    """
    # SIN(VO VA FREQ TD THETA)

    def __init__(self, vo, va, freq, td=0., theta=0., phi=0.):
        self.vo = vo
        self.va = va
        self.freq = freq
        self.td = td
        self.theta = theta
        self.phi = phi
        self._type = "V"

    def __call__(self, time):
        """Evaluate the sine function at the given time."""
        if time is None:
            time = 0
        if time < self.td:
            return self.vo + self.va*math.sin(math.pi*self.phi/180.)
        else:
            return self.vo + self.va * math.exp((self.td - time)*self.theta) \
                   * math.sin(2*math.pi*self.freq*(time - self.td) + \
                              math.pi*self.phi/180.)

    def __str__(self):
        return "type=sin " + \
            self._type.lower() + "o=" + str(self.vo) + " " + \
            self._type.lower() + "a=" + str(self.va) + \
            " freq=" + str(self.freq) + " theta=" + str(self.theta) + \
            " td=" + str(self.td)


class exp(object):
    """Exponential wave

    .. image:: images/elem/exp.svg

    Mathematically, it is described by the equations:

    * :math:`0 \\le t < TD1`:

    .. math::

        f(t) = V1

    * :math:`TD1 < t < TD2`

    .. math::

        f(t) = V1+(V2-V1) \\cdot \\left[1-\\exp
               \\left(-\\frac{t-TD1}{TAU1}\\right)\\right]

    * :math:`t > TD2`

    .. math::

        f(t) = V1+(V2-V1) \\cdot \\left[1-\\exp
               \\left(-\\frac{t-TD1}{TAU1}\\right)\\right]+(V1-V2) \\cdot
               \\left[1-\\exp \\left(-\\frac{t-TD2}{TAU2}\\right)\\right]

    **Parameters:**

    v1 : float
        Initial value.
    v2 : float
        Pulsed value.
    td1 : float
        Rise delay time in seconds.
    td2 : float
        Fall delay time in seconds.
    tau1 : float
        Rise time constant in seconds.
    tau2 : float
        Fall time constant in seconds.
    """
    # EXP(V1 V2 TD1 TAU1 TD2 TAU2)

    def __init__(self, v1, v2, td1, tau1, td2, tau2):
        self.v1 = v1
        self.v2 = v2
        self.td1 = td1
        self.tau1 = tau1
        self.td2 = td2
        self.tau2 = tau2
        self._type = "V"

    def __call__(self, time):
        """Evaluate the exponential function at the given time."""
        if time is None:
            time = 0
        if time < self.td1:
            return self.v1
        elif time < self.td2:
            return self.v1 + (self.v2 - self.v1) * \
                   (1 - math.exp(-1*(time - self.td1)/self.tau1))
        else:
            return self.v1 + (self.v2 - self.v1) * \
                   (1 - math.exp(-1*(time - self.td1)/self.tau1)) + \
                   (self.v1 - self.v2)*(1 - math.exp(-1*(time - self.td2)/self.tau2))

    def __str__(self):
        return "type=exp " + \
            self._type.lower() + "1=" + str(self.v1) + " " + \
            self._type.lower() + "2=" + str(self.v2) + \
            " td1=" + str(self.td1) + " td2=" + str(self.td2) + \
            " tau1=" + str(self.tau1) + " tau2=" + str(self.tau2)


class sffm(object):
    """Single-Frequency FM (SFFM) waveform

    .. image:: images/elem/fm.svg

    Mathematically, it is described by the equations:

    * :math:`0 \\le t \\le t_D`:

    .. math::

        f(t) = V_O

    * :math:`t > t_D`

    .. math::

        f(t) = V_O + V_A \\cdot \\sin \\left[2\\pi f_C (t - t_D) + MDI
               \\sin \\left[2 \\pi f_S (t - t_D) \\right] \\right]

    **Parameters:**

    vo : float
        Offset in Volt or Ampere.
    va : float
        Amplitude in Volt or Ampere.
    fc : float
        Carrier frequency in Hz.
    mdi : float
        Modulation index.
    fs : float
        Signal frequency in HZ.
    td : float
        Time delay before the signal begins, in seconds.
    """
    # SFFM(VO VA FC MDI FS)

    def __init__(self, vo, va, fc, mdi, fs, td):
        self.vo = vo
        self.va = va
        self.fc = fc
        self.mdi = mdi
        self.fs = fs
        self.td = td
        self._type = "V"

    def __call__(self, time):
        """Evaluate the SFFM function at the given time."""
        if time is None:
            time = 0
        if time <= self.td:
            return self.vo
        else:
            return self.vo + self.va*math.sin(2*math.pi*self.fc*(time - self.td) +
                                              self.mdi*math.sin(2*math.pi*self.fs*
                                                                (time - self.td))
                                              )

    def __str__(self):
        return "type=sffm vo=%g va=%g fc=%g mdi=%g fs=%g td=%g" % \
                (self.vo, self.va, self.fc, self.mdi, self.fs, self.td)

class am(object):
    """Amplitude Modulated (AM) waveform

    .. image:: images/elem/am.svg

    Mathematically, it is described by the equations:

    * :math:`0 \\le t \\le t_D`:

    .. math::

        f(t) = O

    * :math:`t > t_D`

    .. math::

        f(t) = SA \\cdot \\left[OC + \\sin \\left[2\\pi f_m (t - t_D) \\right]
               \\right] \\cdot \\sin \\left[2 \\pi f_c (t - t_D) \\right]

    **Parameters:**

    sa : float
        Signal amplitude in Volt or Ampere.
    fc : float
        Carrier frequency in Hertz.
    fm : float
        Modulation frequency in Hertz.
    oc : float
        Offset constant, setting the absolute magnitude of the modulation.
    td : float
        Time delay before the signal begins, in seconds.
    """
    # AM(sa oc fm fc <td>)

    def __init__(self, sa, fc, fm, oc, td):
        self.sa = sa
        self.fc = fc
        self.fm = fm
        self.oc = oc
        self.td = td
        self._type = "V"

    def __call__(self, time):
        """Evaluate the AM function at the given time."""
        if time is None:
            time = 0
        if time <= self.td:
            return 0.
        else:
            return self.sa*(self.oc + math.sin(2*math.pi*self.fm*
                                               (time - self.td)))* \
                   math.sin(2*math.pi*self.fc*(time - self.td))

    def __str__(self):
        return "type=am sa=%g oc=%g fm=%g fc=%g td=%g" % \
                (self.sa, self.oc, self.fm, self.fc, self.td)

class pwl(object):
    """Piece-Wise Linear (PWL) waveform

    .. image:: images/elem/pwl.svg

    A piece-wise linear waveform is defined by a sequence of points
    :math:`(x_i, y_i)`.

    Please supply the abscissa values :math:`\\{x\\}_i` in the vector
    ``x``, the ordinate values :math:`\\{y\\}_i` in the vector ``y``,
    separately.


    **Parameters:**

    x : sequence-like
        The abscissa values of the interpolation points.
    y : sequence-like
        The ordinate values of the interpolation points.
    repeat : boolean, optional
        Whether the waveform should be repeated after its end. If set to
        ``True``, ``repeat_time`` also needs to be set to define when the
        repetition begins. Defaults to ``False``.
    repeat_time : float, optional
        In case the waveform is set to be repeated, setting the ``repeat`` flag
        above, the parameter, defined in seconds, set the first time instant at
        which the waveform repetition happens.
    td : float, optional
        Time delay before the signal begins, in seconds. Defaults to zero.

    **Example:**

    The following code::

        import ahkab
        import numpy as np
        import pylab as plt
        # vs = (x1, y1, x2, y2, x3, y3 ...)
        vs = (60e-9, 0, 120e-9, 0, 130e-9, 5, 170e-9, 5, 180e-9, 0)
        x, y = vs[::2], vs[1::2]
        fun = ahkab.time_functions.pwl(x, y, repeat=1, repeat_time=60e-9, td=0)
        myg = np.frompyfunc(fun, 1, 1)
        t = np.linspace(0, 5e-7, 2000)
        plt.plot(t, myg(t), lw=3)
        plt.xlabel('Time [s]'); plt.ylabel('Arbitrary units []')


    Produces:

    .. plot::

        import ahkab
        import numpy as np
        import pylab as plt
        vs = (60e-9, 0, 120e-9, 0, 130e-9, 5, 170e-9, 5, 180e-9, 0)
        x, y = vs[::2], vs[1::2]
        fun = ahkab.time_functions.pwl(x, y, repeat=1, repeat_time=60e-9, td=0)
        myg = np.frompyfunc(fun, 1, 1)
        t = np.linspace(0, 5e-7, 2000)
        plt.figure(figsize=(6, 3)); plt.grid()
        plt.plot(t, myg(t), lw=3)
        plt.xlabel('Time [s]'); plt.ylabel('Arbitrary units []')
        plt.tight_layout()


    """

    def __init__(self, x, y, repeat=False, repeat_time=0, td=0):
        self.x = x
        self.y = y
        self.repeat = repeat
        self.repeat_time = repeat_time
        if self.repeat_time == max(x):
            self.repeat_time = 0
        self.td = td
        self._type = "V"
        self._f = InterpolatedUnivariateSpline(self.x, self.y, k=1)

    def __call__(self, time):
        """Evaluate the PWL function at the given time."""
        time = self._normalize_time(time)
        return self._f(time)

    def _normalize_time(self, time):
        if time is None:
            time = 0
        if time <= self.td:
            time = 0
        elif time > self.td:
            time = time - self.td
            if self.repeat:
                if time > max(self.x):
                    time = (time - max(self.x)) % \
                           (max(self.x) - self.repeat_time) + \
                           self.repeat_time
                else:
                    pass
        return time

    def __str__(self):
        pwl_str = "type=pwl"
        tv = " "
        for x, y in zip(self.x, self.y):
            tv += "%g %g "
        pwl_str += tv
        if self.td:
            pwl_str += "td=%g " % self.td
        if self.repeat:
            pwl_str = "RPT=%g " % self.repeat_time
        return pwl_str[:-1]

