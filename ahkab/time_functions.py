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
are intefacing to the simulator through Python. Please see the dedicated section
below.

Classes defined in this module
------------------------------

.. autosummary::
    pulse
    sin
    exp
    sffm
    am

Defining custom time functions
------------------------------

Defining a custom time function is easy, all you need is:

* An object with a ``value(self, time)`` method.

The simulator will call ``value(self, time)`` of the class instance you provide
at every time step in time-based simulations. It expects to receive as return
value a ``float``, corresponding to the value of the voltage applied by the
voltage source, in Volt, if the custom time function was passed to
:class:`VSource`, or to the value of the current flowing through the current
source, if the custom time function was passed to :class:`ISource`.

The standard notation applies.

Module reference
----------------
"""

from __future__ import (unicode_literals, absolute_import,
                        division, print_function)

import math

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

class pulse:
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

    def value(self, time):
        """Evaluate the pulse function at the given time."""
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


class sin:
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

    def value(self, time):
        """Evaluate the sine function at the given time."""
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


class exp:
    """Exponential time function

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

    def value(self, time):
        """Evaluate the exponential function at the given time."""
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


class sffm:
    """Single-Frequency FM time function

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

    def value(self, time):
        """Evaluate the SFFM function at the given time."""
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

class am:
    """AM time function

    .. image:: images/elem/am.svg

    Mathematically, it is described by the equations:

    * :math:`0 \\le t \\le t_D`:

    .. math::

        f(t) = O

    * :math:`t > t_D`

    .. math::

        f(t) = SA \\cdot \\left[OC + \\sin \\left[2\\pi f_m (t - t_D) \\right]
               \\cdot \\sin \\left[2 \\pi f_c (t - t_D) \\right]

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

    def value(self, time):
        """Evaluate the AM function at the given time."""
        if time <= self.td:
            return 0.
        else:
            return self.sa*(self.oc + math.sin(2*math.pi*self.fm*
                                               (time - self.td)))* \
                   math.sin(2*math.pi*self.fc*(time - self.td))

    def __str__(self):
        return "type=am sa=%g oc=%g fm=%g fc=%g td=%g" % \
                (self.sa, self.oc, self.fm, self.fc, self.td)

