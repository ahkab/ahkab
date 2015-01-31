# -*- coding: iso-8859-1 -*-
# devices.py
# Devices for simulation
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
This module contains several basic element classes and time functions.

Introduction
------------

While they may be instantiated directly by the user, notice that the
main ``ahkab`` module provides convenience functions to instantiate
and connect into a circuit instance all of the following devices.

Of more interest to the end user are the time function classes,
which the user will have to instantiate to provide a time-varying
characteristic to independent sources.

Notice that both time functions and circuit elements are not restricted
to those provided here, the user is welcome to provide his own.

While implementing a new component requires some understanding of the internals
of the simultor and it is expected to be less common, implementing a custom time
function is easy and common practice, as long as you are intefacing to the
simulator through Python.

Both of the two cases have a dedicated section below.

Classes defined in this module
------------------------------

Elements
========

.. autosummary::
    ISource
    VSource
    Resistor
    Capacitor
    Inductor
    InductorCoupling
    EVSource
    GISource
    HVSource
    FISource

Time functions
===============

.. autosummary::
    pulse
    sin
    exp

Defining new elements and subclassing ``Component``
---------------------------------------------------

We recommend to subclass :class:`ahkab.devices.Component` if you intend to
define a new element.

The general form of a (possibly nonlinear) element class is described in the
following.

Required attributes and methods
===============================

The class must provide:

1. Element terminals:

::

    elem.n1 # the anode of the element
    elem.n2 # the cathode of the element

.. note:: a positive current is a current that flows into the anode and out of
    the cathode. This convention is used throughout the simulator.

2. ``elem.get_ports()``

This method must return a tuple of pairs of nodes.

Eg.

::

    ((na, nb), (nc, nd), (ne, nf), ... )

Each pair of nodes is used to determine a voltage that has effect on the
current.

For example, the source-referred model of an nmos may provide:

::

    ((n_gate, n_source), (n_drain, n_source))

The positive terminal is the first.

From that, the calling method builds a voltage vector corresponding to the
ports vector:

::

    voltages_vector = ( Va-Vb, Vc-Vd, Ve-Vf, ...)

That's passed to:

3. ``elem.i(voltages_vector, time)``

It returns the current flowing into the element if the voltages specified in
the voltages_vector are applied to its ports, at the time given.

4. ``elem.g(voltages_vector, port_index, time)``

similarly returns the differential transconductance between the port at
position ``port_index`` in the ``ports_vector`` (see point **2** above)
and the element output current, when the operating point is specified by
the voltages in the ``voltages_vector``.

5. ``elem.is_nonlinear``

A non linear element must have a ``elem.is_nonlinear`` field set to True.

6. ``elem.is_symbolic``

This boolean flag is used to know whether the element should be treated
symbolically by the ymbolic solver or not. It is meant to be toggled
by the user at will.

7. Every element should have a ``get_netlist_elem_line(self, nodes_dict)``
allowing the element to print a netlist entry that parses to itself.

Recommended attributes and methods
==================================

1. A non linear element may have a list/tuple of the same length of its
``ports_vector`` in which there are the recommended guesses for DC and OP
analyses.

Eg. ``Vgs`` is set to ``Vt0`` in mosfets.

This is obviously useless for linear devices.

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

import numpy as np
import math

from . import constants
from . import printing

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
               },
               )
}}


class Component(object):

    """Base Component class.

    This component is not meant for direct use, rather all other (simple)
    components are a subclass of this element.

    """

    def __init__(self, part_id=None, n1=None, n2=None, is_nonlinear=False, is_symbolic=True, value=None):
        self.part_id = part_id
        self.n1 = n1
        self.n2 = n2
        self.value = value
        self.is_nonlinear = is_nonlinear
        self.is_symbolic = is_symbolic

    #   Used by `get_netlist_elem_line` for value
    def __str__(self):
        return str(self.value)

    #   must be called to define the element!
    def set_char(self, i_function=None, g_function=None):
        if i_function:
            self.i = i_function
        if g_function:
            self.g = g_function

    def g(self, v):
        return 1./self.value

    def i(self, v):
        return 0

    def get_netlist_elem_line(self, nodes_dict):
        """A netlist line that, parsed, evaluates to the same instance

        **Parameters:**

        nodes_dict : dict
            The nodes dictionary of the circuit, so that the method
            can convert its internal node IDs to the corresponding
            external ones.

        **Returns:**

        ntlst_line : string
            The netlist line.
        """
        return "%s %s %s %g" % (self.part_id, nodes_dict[self.n1],
                                nodes_dict[self.n2], self.value)


class Resistor(Component):
    """A resistor.

    .. image:: images/elem/resistor.svg

    **Parameters:**

    part_id : string
        The unique identifier of this element. The first letter should be
        ``'R'``.
    n1 : int
        *Internal* node to be connected to the anode.
    n2 : int
        *Internal* node to be connected to the cathode.
    value : float
        Resistance in ohms.

     """
    #
    #             /\    /\    /\
    #     n1 o---+  \  /  \  /  \  +---o n2
    #                \/    \/    \/
    #
    def __init__(self, part_id, n1, n2, value):
        self.part_id = part_id
        self._value = value
        self._g = 1./value
        self.is_nonlinear = False
        self.is_symbolic = True
        self.n1 = n1
        self.n2 = n2

    @property
    def g(self, v=0, time=0):
        return self._g

    @g.setter
    def g(self, g):
        self._g = g
        self._value = 1./g

    @property
    def value(self, v=0, time=0):
        return self._value

    @value.setter
    def value(self, value):
        self._value = value
        self._g = 1./value

    def i(self, v, time=0):
        return 0

    def get_op_info(self, ports_v):
        vn1n2 = float(ports_v[0][0])
        in1n2 = float(ports_v[0][0] / self.value)
        power = float(ports_v[0][0] ** 2 / self.value)
        arr = [
            [self.part_id.upper(), "V(n1-n2):", vn1n2, "[V]", "I(n2-n1):", in1n2, "[A]", "P:", power, "[W]"]]
        strarr = printing.table_setup(arr)
        return strarr

    def print_op_info(self, ports_v):
        print(self.get_op_info(ports_v))


class Capacitor(Component):
    """A capacitor.

    .. image:: images/elem/capacitor.svg

    **Parameters:**

    part_id : string
        The unique identifier of this element. The first letter should be
        ``'C'``.
    n1 : int
        *Internal* node to be connected to the anode.
    n2 : int
        *Internal* node to be connected to the cathode.
    value : float
        The capacitance in Farads.
    ic : float
        The initial condition (IC) to be used for time-based simulations,
        such as TRAN analyses, when requested, expressed in Volt.

    """
    #
    #               |  |
    #               |  |
    #     n1 o------+  +-------o n2
    #               |  |
    #               |  |
    #
    def __init__(self, part_id, n1, n2, value, ic=None):
        self.part_id = part_id
        self.value = value
        self.n1 = n1
        self.n2 = n2
        self.ic = ic
        self.is_nonlinear = False
        self.is_symbolic = True

    def g(self, v, time=0):
        return 0

    def i(self, v, time=0):
        return 0

    def d(self, v, time=0):
        return self.value

    def get_op_info(self, ports_v):
        vn1n2 = float(ports_v[0][0])
        qn1n2 = float(ports_v[0][0] * self.value)
        energy = float(.5 * ports_v[0][0] ** 2 * self.value)
        arr = [
            [self.part_id.upper(), "V(n1-n2):", vn1n2, "[V]", "Q:", qn1n2, "[C]", "E:", energy, "[J]"]]
        strarr = printing.table_setup(arr)
        return strarr

    def print_op_info(self, ports_v):
        print(self.get_op_info(ports_v))


class Inductor(Component):
    """An inductor.

    .. image:: images/elem/inductor.svg

    **Parameters:**

    part_id : string
        The unique identifier of this element. The first letter should be
        ``'L'``.
    n1 : int
        *Internal* node to be connected to the anode.
    n2 : int
        *Internal* node to be connected to the cathode.
    value : float
        The inductance in Henry.
    ic : float
        The initial condition (IC) to be used for time-based simulations,
        such as TRAN analyses, when requested, expressed in Ampere.

    """
    #
    #             L
    #  n1 o----((((((((----o n2
    #
    #
    def __init__(self, part_id, n1, n2, value, ic=None):
        self.value = value
        self.n1 = n1
        self.n2 = n2
        self.part_id = part_id
        self.ic = ic
        self.coupling_devices = []
        self.is_nonlinear = False
        self.is_symbolic = True


class InductorCoupling(Component):
    """Coupling between two inductors.

    .. image:: images/elem/mutual_inductors.svg

    This element is used to simulate the coupling between two inductors,
    such as in the case of a transformer.

    Notice that turn ratio and the inductance ratio are linked by the
    relationship:

    .. math::

        \\frac{L_1}{L_2} = \\left(\\frac{N_1}{N_2}\\right)^2

    **Parameters:**

    part_id : string
        The unique identifier of this element. The first letter should be
        ``'K'``.
    L1 : string
        The ``part_id`` of the first inductor to be coupled.
    L2 : string
        The ``part_id`` of the second inductor to be coupled.
    K : float
        The coupling coefficient between the two windings.
    M : float
        The mutual inductance between the windings, it is equal to
        :math:`K\\sqrt(L_1L2)`, where :math:`L_1` and :math:`L_2` are the
        values of the two inductors ``L1`` and ``L2``.
    """
    # K1 L1 L2 k=<float>
    # M = sqrt(L1elem.value * L2elem.value) * Kvalue
    def __init__(self, part_id, L1, L2, K, M):
        self.part_id = part_id
        self.L1 = L1
        self.L2 = L2
        self.M = M
        self.K = K
        self.is_nonlinear = False
        self.is_symbolic = True

    def __str__(self):
        return "%s %s %g" % (self.L1, self.L2, self.value)

    def get_other_inductor(self, Lselected):
        Lret = None
        if Lselected.upper() == self.L1.upper():
            Lret = self.L2
        elif Lselected.upper() == self.L2.upper():
            Lret = self.L1
        if Lret is None:
            raise Exception("Mutual inductors bug.")
        return Lret

    def get_netlist_elem_line(self, nodes_dict):
        """A netlist line that, parsed, evaluates to the same instance

        **Parameters:**

        nodes_dict : dict
            The nodes dictionary of the circuit, so that the method
            can convert its internal node IDs to the corresponding
            external ones.

        **Returns:**

        ntlst_line : string
            The netlist line.
        """
        return "%s %s %s %g" % (self.part_id, self.L1, self.L2, self.K)


#
# SOURCES
#


class ISource(Component):
    """An ideal current source.

    .. image:: images/elem/isource.svg

    Defaults to a DC current source.

    To implement a time-varying source:

    * set ``_time_function`` to an appropriate instance having a
      ``value(self, time)`` method,
    * set ``is_timedependent`` to ``True``.

    **Parameters:**

    part_id : string
        The unique identifier of this element. The first letter should be
        ``'I'``.
    n1 : int
        *Internal* node to be connected to the anode.
    n2 : int
        *Internal* node to be connected to the cathode.
    dc_value : float
        DC voltage in Ampere.
    ac_value : complex float, optional
        AC current in Ampere. Defaults to no AC characteristics,
        ie :math:`I(\\omega) = 0 \\;\\;\\forall \\omega > 0`.

    """
    def __init__(self, part_id, n1, n2, dc_value=None, ac_value=0):
        self.part_id = part_id
        self.dc_value = dc_value
        self.abs_ac = np.abs(ac_value) if ac_value else None
        self.arg_ac = np.angle(ac_value) if ac_value else None
        self.n1 = n1
        self.n2 = n2
        self.is_nonlinear = False
        self.is_symbolic = True
        self.is_timedependent = False
        self._time_function = None

    def __str__(self):
        rep = ""
        if self.dc_value is not None:
            rep = rep + "type=idc value=" + str(self.dc_value) + " "
        if self.abs_ac is not None:
            rep = rep + "iac=" + \
                str(self.abs_ac) + " " + "arg=" + str(self.arg_ac) + " "
        if self.is_timedependent:
            rep = rep + str(self._time_function)
        return rep

    def I(self, time=None):
        """Evaluate the current forced by the current source.

        If ``time`` is not supplied, or if it is set to ``None``, or if the
        source is only specified for DC, returns ``dc_value``.

        **Parameters:**

        time : float or None, optional
            The time at which the current is evaluated, if any.

        **Returns:**

        I : float
            The current, in Ampere.

        .. note::

            This simulator uses passive convention:
            A positive currents flows in a element into the positive node and
            out of the negative node
        """
        if not self.is_timedependent or (self._time_function == None) or \
            (time == None and self.dc_value is not None):
            return self.dc_value
        else:
            return self._time_function.value(time)

    def get_netlist_elem_line(self, nodes_dict):
        """A netlist line that, parsed, evaluates to the same instance

        **Parameters:**

        nodes_dict : dict
            The nodes dictionary of the circuit, so that the method
            can convert its internal node IDs to the corresponding
            external ones.

        **Returns:**

        ntlst_line : string
            The netlist line.
        """
        rep = ""
        rep += "%s %s %s " % (self.part_id, nodes_dict[self.n1],
                             nodes_dict[self.n2])
        if self.dc_value is not None:
            rep = rep + "type=idc value=" + str(self.dc_value) + " "
        if self.abs_ac is not None:
            rep = rep + "iac=" + \
                str(self.abs_ac) + " " + "arg=" + str(self.arg_ac) + " "
        if self.is_timedependent:
            rep = rep + str(self._time_function)
        return rep


class VSource(Component):
    """An ideal voltage source.

    .. image:: images/elem/vsource.svg

    Defaults to a DC voltage source.

    To implement a time-varying source:

    * set ``_time_function`` to an appropriate instance having a
      ``value(self, time)`` method,
    * set ``is_timedependent`` to ``True``.

    **Parameters:**

    part_id : string
        The unique identifier of this element. The first letter should be
        ``'V'``.
    n1 : int
        *Internal* node to be connected to the anode.
    n2 : int
        *Internal* node to be connected to the cathode.
    dc_value : float
        DC voltage in Volt.
    ac_value : complex float, optional
        AC voltage in Volt. Defaults to no AC characteristics,
        ie :math:`V(\\omega) = 0 \\;\\;\\forall \\omega > 0`.

    """

    def __init__(self, part_id, n1, n2, dc_value, ac_value=0):
        self.part_id = part_id
        self.dc_value = dc_value
        self.n1 = n1
        self.n2 = n2
        self.abs_ac = np.abs(ac_value) if ac_value else None
        self.arg_ac = np.angle(ac_value) if ac_value else None
        self.is_nonlinear = False
        self.is_symbolic = True
        self.is_timedependent = False
        self._time_function = None
        if dc_value is not None:
            self.dc_guess = [self.dc_value]

    def __str__(self):
        rep = ""
        if self.dc_value is not None:
            rep = rep + "type=vdc value=" + str(self.dc_value) + " "
        if self.abs_ac is not None:
            #   TODO:   netlist parser doesn't accept `arg=` from `self.arg_ac`
            rep = rep + "vac=" + str(self.abs_ac) + " "
        if self.is_timedependent:
            rep = rep + str(self._time_function)
        return rep

    def V(self, time=None):
        """Evaluate the voltage applied by the voltage source.

        If ``time`` is not supplied, or if it is set to ``None``, or if the
        source is only specified for DC, returns ``dc_value``.

        **Parameters:**

        time : float or None, optional
            The time at which the voltage is evaluated, if any.

        **Returns:**

        V : float
            The voltage, in Volt.
        """

        if not self.is_timedependent or \
            (self._time_function is None) or \
                (time is None and self.dc_value is not None):
            return self.dc_value
        else:
            return self._time_function.value(time)

    def get_netlist_elem_line(self, nodes_dict):
        """A netlist line that, parsed, evaluates to the same instance

        **Parameters:**

        nodes_dict : dict
            The nodes dictionary of the circuit, so that the method
            can convert its internal node IDs to the corresponding
            external ones.

        **Returns:**

        ntlst_line : string
            The netlist line.
        """
        rep = ""
        rep += "%s %s %s " % (self.part_id, nodes_dict[self.n1],
                             nodes_dict[self.n2])
        if self.dc_value is not None:
            rep = rep + "type=vdc value=" + str(self.dc_value) + " "
        if self.abs_ac is not None:
            #   TODO:   netlist parser doesn't accept `arg=` from `self.arg_ac`
            rep = rep + "vac=" + str(self.abs_ac) + " "
        if self.is_timedependent:
            rep = rep + str(self._time_function)
        return rep


class EVSource(Component):
    """Linear voltage-controlled voltage source

    .. image:: images/elem/vcvs.svg

    The source port is an open circuit, the destination port is an ideal
    voltage source.

    Mathematically, it is equivalent to the following:

    .. math::

        \\left\\{
        \\begin{array}{ll}
            I_s = 0\\\\
            Vn_1 - Vn_2 = \\alpha * (Vsn_1 - Vsn_2)
        \\end{array}
        \\right.

    Where :math:`I_s` is the current at the source port and the remaining
    symbols are explained in the Parameters section below.

    **Parameters:**

    n1 : int
        *Internal* node to be connected to the anode of the output port.
    n2 : int
        *Internal* node to be connected to the cathode of the output port.
    value : float
        Proportionality constant :math:`\\alpha` between the voltages.
    sn1 : int
        *Internal* node to be connected to the anode of the source (sensing)
        port.
    sn2 : int
        *Internal* node to be connected to the cathode of the source
        (sensing) port.
    """

    def __init__(self, part_id, n1, n2, value, sn1, sn2):
        self.part_id = part_id
        self.n1 = n1
        self.n2 = n2
        self.alpha = value
        self.sn1 = sn1
        self.sn2 = sn2
        self.is_nonlinear = False
        self.is_symbolic = True

    def __str__(self):
        return "alpha=%s" % self.alpha

    def get_netlist_elem_line(self, nodes_dict):
        """A netlist line that, parsed, evaluates to the same instance

        **Parameters:**

        nodes_dict : dict
            The nodes dictionary of the circuit, so that the method
            can convert its internal node IDs to the corresponding
            external ones.

        **Returns:**

        ntlst_line : string
            The netlist line.
        """
        return "%s %s %s %s %s %g" % (self.part_id, nodes_dict[self.n1],
                                nodes_dict[self.n2], nodes_dict[self.sn1],
                                nodes_dict[self.sn2], self.alpha)


class GISource(Component):

    """Linear voltage controlled current source

    .. image:: images/elem/vccs.svg

    The source port is an open circuit, the output port is an ideal current
    source:

    .. math::

        \\left\\{
        \\begin{array}{ll}
            I_s = 0\\\\
            I_o = \\alpha \\cdot (V(sn_1) - V(sn_2))
        \\end{array}
        \\right.


    Where :math:`I_s` is the current at the source port and :math:`I_o` is the
    current at the output port.
    The remaining symbols are explained in the Parameters section below.

    .. note::

        This simulator uses the passive convention: a positive current flows
        into the element through the anode and exits through the cathode.

    **Parameters:**

    n1 : int
        *Internal* node to be connected to the anode of the output port.
    n2 : int
        *Internal* node to be connected to the cathode of the output port.
    value : float
        Proportionality constant :math:`\\alpha` between the sense voltage and
        the output current, in Ampere/Volt.
    sn1 : int
        *Internal* node to be connected to the anode of the source (sensing)
        port.
    sn2 : int
        *Internal* node to be connected to the cathode of the source
        (sensing) port.

    """

    def __init__(self, part_id, n1, n2, value, sn1, sn2):
        self.part_id = part_id
        self.n1 = n1
        self.n2 = n2
        self.alpha = value
        self.sn1 = sn1
        self.sn2 = sn2
        self.is_nonlinear = False
        self.is_symbolic = True

    def __str__(self):
        return "value=%s" % self.alpha

    def get_netlist_elem_line(self, nodes_dict):
        """A netlist line that, parsed, evaluates to the same instance

        **Parameters:**

        nodes_dict : dict
            The nodes dictionary of the circuit, so that the method
            can convert its internal node IDs to the corresponding
            external ones.

        **Returns:**

        ntlst_line : string
            The netlist line.
        """
        return "%s %s %s %s %s %g" % (self.part_id, nodes_dict[self.n1],
                                nodes_dict[self.n2], nodes_dict[self.sn1],
                                nodes_dict[self.sn2], self.alpha)


class HVSource(Component):  # H

    """Linear current-controlled voltage source

    .. image:: images/elem/ccvs.svg

    The source port is an existing voltage source, used to sense the current
    controlling the voltage source connected to the destination port.

    Mathematically, it is equivalent to the following:

    .. math::

        \\left\\{
        \\begin{array}{ll}
            V(sn_1) - V(sn_2) = V_S \\\\
            Vn_1 - Vn_2 = \\alpha * I[V_s]
        \\end{array}
        \\right.

    Where :math:`I[V_s]` is the current flowing in the source port, :math:`V_s`
    is the voltage applied between the nodes :math:`sn_1` and :math:`sn_2`.
    The remaining symbols are explained in the Parameters section below.

    .. note::

        This simulator uses the passive convention: a positive current flows
        into the element through the anode and exits through the cathode.

    **Parameters:**

    n1 : int
        *Internal* node to be connected to the anode of the output port.
    n2 : int
        *Internal* node to be connected to the cathode of the output port.
    value : float
        Proportionality constant :math:`\\alpha` between the sense current and
        the output voltage, in V/A.
    source_id : string
        ``part_id`` of the current-sensing voltage source, eg. ``'V1'`` or
        ``'VSENSE'``.

    """

    def __init__(self, part_id, n1, n2, value, source_id):
        self.part_id = part_id
        self.n1 = n1
        self.n2 = n2
        self.alpha = value
        self.source_id = source_id
        self.is_nonlinear = False
        self.is_symbolic = True

    def get_netlist_elem_line(self, nodes_dict):
        """A netlist line that, parsed, evaluates to the same instance

        **Parameters:**

        nodes_dict : dict
            The nodes dictionary of the circuit, so that the method
            can convert its internal node IDs to the corresponding
            external ones.

        **Returns:**

        ntlst_line : string
            The netlist line.
        """
        return "%s %s %s %s %g" % (self.part_id, nodes_dict[self.n1],
                                nodes_dict[self.n2], self.source_id,
                                self.alpha)


class FISource(Component):
    """Linear current-controlled current source

    .. image:: images/elem/cccs.svg

    This element implements a current source whose current value is controlled
    by the current flowing in a current source, which acts as the "sensing"
    element.

    Mathematically:

    .. math::


        \\left\\{
        \\begin{array}{ll}
            V(sn_1) - V(sn_2) = V_S \\\\
            I_o = \\alpha \\cdot I_s
        \\end{array}
        \\right.


    Where :math:`V_s` is the voltage forced at the source port by the sensing
    element and :math:`I_o` is the current at the output port.
    The remaining symbols are explained in the Parameters section below.

    .. note::

        This simulator uses the passive convention: a positive current flows
        into the element through the anode and exits through the cathode.

    **Parameters:**

    n1 : int
        *Internal* node to be connected to the anode of the output port.
    n2 : int
        *Internal* node to be connected to the cathode of the output port.
    value : float
        Proportionality constant :math:`\\alpha` between the sense current and
        the output current.
    source_id : string
        ``part_id`` of the sensing voltage source, eg. ``'V1'`` or ``'VSENSE'``.

    """
    #F
    def __init__(self, part_id, n1, n2, value, source_id):
        self.part_id = part_id
        self.n1 = n1
        self.n2 = n2
        self.source_id = source_id
        self.alpha = value
        self.is_nonlinear = False
        self.is_symbolic = True

    def get_netlist_elem_line(self, nodes_dict):
        """A netlist line that, parsed, evaluates to the same instance

        **Parameters:**

        nodes_dict : dict
            The nodes dictionary of the circuit, so that the method
            can convert its internal node IDs to the corresponding
            external ones.

        **Returns:**

        ntlst_line : string
            The netlist line.
        """
        return "%s %s %s %s %g" % (self.part_id, nodes_dict[self.n1],
                                nodes_dict[self.n2], self.source_id,
                                self.alpha)


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

        f(t) = v_o + v_a \\exp\\left[-(t - t_d)\,\\mathrm{THETA}\\right] \\sin\\left[2 \\pi f (t - t_d) + \\pi \\phi/180\\right]

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
