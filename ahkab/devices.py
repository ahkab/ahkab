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
This module contains several basic element classes.

Introduction
------------

While they may be instantiated directly by the user, notice that the
main ``ahkab`` module provides convenience functions to instantiate
and connect into a circuit instance all of the following devices.

Notice that the circuit elements are not restricted to those provided here, the
user is welcome to provide his own. Please see the dedicated section below.

Classes defined in this module
------------------------------

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

Module reference
----------------

"""

from __future__ import (unicode_literals, absolute_import,
                        division, print_function)

import numpy as np

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
        """Information regarding the Operating Point (OP)

        **Parameters:**

        ports_v : list of lists
            The parameter is to be set to ``[[v]]``, where ``v`` is the voltage
            applied to the resistor terminals.

        **Returns:**

        op_keys : list of strings
            The labels corresponding to the numeric values in ``op_info``.
        op_info : list of floats
            The values corresponding to ``op_keys``.
        """
        vn1n2 = float(ports_v[0][0])
        in1n2 = float(ports_v[0][0]/self.value)
        power = float(ports_v[0][0] ** 2 / self.value)
        op_keys = ['Part ID', u"R [\u2126]", "V(n1,n2) [V]", "I(n1->n2) [A]", "P [W]"]
        op_info = [self.part_id.upper(), self.value, vn1n2, in1n2, power]
        return op_keys, op_info

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
        """Information regarding the Operating Point (OP)

        **Parameters:**

        ports_v : list of lists
            The parameter is to be set to ``[[v]]``, where ``v`` is the voltage
            applied to the capacitor terminals.

        **Returns:**

        op_keys : list of strings
            The labels corresponding to the numeric values in ``op_info``.
        op_info : list of floats
            The values corresponding to ``op_keys``.
        """
        vn1n2 = float(ports_v[0][0])
        qn1n2 = float(ports_v[0][0] * self.value)
        energy = float(.5 * ports_v[0][0] ** 2 * self.value)
        op_keys = ['Part ID', "V(n1-n2) [V]", "Q [C]", "E [J]"]
        op_info = [self.part_id.upper(), vn1n2, qn1n2, energy]
        return op_keys, op_info

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

    def get_op_info(self, ports_v, current):
        """Information regarding the Operating Point (OP)

        **Parameters:**

        ports_v : list of lists
            The parameter is to be set to ``[[v]]``, where ``v`` is the voltage
            applied to the inductor terminals.
        current : float
            The current flowing in the inductor, positive currents flow in ``n1``
            and out of ``n2``.

        **Returns:**

        op_keys : list of strings
            The labels corresponding to the numeric values in ``op_info``.
        op_info : list of floats
            The values corresponding to ``op_keys``.
        """
        vn1n2 = float(ports_v[0][0])
        energy = .5 * self.value * current**2
        op_keys = ['Part ID', u"\u03d5(n1,n2) [Wb]", "I(n1->n2) [A]", "E [J]"]
        op_info = [self.part_id.upper(), self.value*current, current, energy]
        return op_keys, op_info


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
        if not self.is_timedependent or (self._time_function is None) or \
            (time is None and self.dc_value is not None):
            return self.dc_value
        else:
            return self._time_function(time)

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

    def get_op_info(self, ports_v):
        """Information regarding the Operating Point (OP)

        **Parameters:**

        ports_v : list of lists
            The parameter is to be set to ``[[v]]``, where ``v`` is the voltage
            applied to the current source terminals.

        **Returns:**

        op_keys : list of strings
            The labels corresponding to the numeric values in ``op_info``.
        op_info : list of floats
            The values corresponding to ``op_keys``.
        """
        vn1n2 = float(ports_v[0][0])
        power = float(ports_v[0][0] * self.I())
        op_keys = ['Part ID', "V(n1-n2) [V]", "I [A]", "P [W]"]
        op_info = [self.part_id.upper(), vn1n2, self.I(), power]
        return op_keys, op_info

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

        if (not self.is_timedependent or\
            self._time_function is None) or \
                (time is None and self.dc_value is not None):
            return self.dc_value
        else:
            return self._time_function(time)

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

    def get_op_info(self, ports_v, current):
        """Information regarding the Operating Point (OP)

        **Parameters:**

        ports_v : list of lists
            The parameter is to be set to ``[[v]]``, where ``v`` is the voltage
            applied to the source terminals.
        current : float
            The current flowing in the voltage source, positive currents flow in
            ``n1`` and out of ``n2``.

        **Returns:**

        op_keys : list of strings
            The labels corresponding to the numeric values in ``op_info``.
        op_info : list of floats
            The values corresponding to ``op_keys``.
        """
        vn1n2 = float(ports_v[0][0])
        power = self.V() * current
        op_keys = ['Part ID', "V(n1,n2) [V]", "I(n1->n2) [A]", "P [W]"]
        op_info = [self.part_id.upper(), self.V(), current, power]
        return op_keys, op_info


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

