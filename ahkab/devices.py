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

__version__ = "0.091"

import numpy
import math
import constants
import printing

"""
Contains various element classes.

General form of a _non linear_ element class:

The class must provide:

1. Element terminals:
elem.n1 # the anode of the element
elem.n2 # the cathode of the element

Notice: a positive current is a current that flows into the anode and out of
the cathode. This convention is used throughout the simulator.

2. elem.get_ports()
This method must return a tuple of pairs of nodes. Something like:
((na, nb), (nc, nd), (ne, nf), ... )

Each pair of nodes is used to determine a voltage that has effect on the
current. I referred to them as a 'port' though it may not be a good idea.

For example, an nmos has:
((n_gate, n_source), (n_drain, n_source))

The positive terminal is the first.

From that, the calling method builds a voltage vector corresponding to the
ports vector:
voltages_vector = ( Va-Vb, Vc-Vd, Ve-Vf, ...)

That's passed to:
3. elem.i(voltages_vector, time)

It returns the current flowing into the element if the voltages specified in
the voltages_vector are applied to its ports, at the time given.

4. elem.g(voltages_vector, port_index, time) is similar, but returns the
differential transconductance between the port at position port_index in the
ports_vector (see 2) and the current, when the operating point is specified by
the voltages in the voltages_vector.

5. elem.is_nonlinear
A non linear element must have a elem.is_nonlinear field set to True.

6. elem.is_symbolic
This boolean flag is used to know whether the element should be treated
symbolically by the ymbolic solver or not. It is meant to be toggled
by the user at will.

Recommended:
1. A non linear element may have a list/tuple of the same length of its
ports_vector in which there are the recommended guesses for dc analysis.
Eg Vgs is set to Vt in mosfets.
This is obviously useless for linear devices.

2. Every element should have a meaningful __str__ method.
It must return a line of paramaters without n1 n2, because a element cannot
know the external names of its nodes. It is used to print the parsed netlist.

"""


class Component(object):

    """Base Component class, also for debugging"""

    def __init__(self, part_id=None, n1=None, n2=None, is_nonlinear=False, is_symbolic=True, value=None):
        self.part_id = part_id
        self.n1 = n1
        self.n2 = n2
        self.value = value
        self.is_nonlinear = is_nonlinear
        self.is_symbolic = is_symbolic

    #   Used by `print_netlist_elem_line` for value
    def __str__(self):
        return str(self.value)

    #   must be called to define the element!
    def set_char(self, i_function=None, g_function=None):
        if i_function:
            self.i = i_function
        if g_function:
            self.g = g_function

    def g(self, v):
        return 1 / self.value

    def i(self, v):
        return 0


class Resistor(Component):
    is_nonlinear = False
    is_symbolic = True

    def __init__(self, part_id='R', n1=None, n2=None, value=None):
        self.part_id = part_id
        self.value = value
        self.n1 = n1
        self.n2 = n2

    def g(self, v, time=0):
        return 1.0 / self.value

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
        print self.get_op_info(ports_v)


class Capacitor(Component):
    is_nonlinear = False
    is_symbolic = True

    def __init__(self, part_id='C', n1=None, n2=None, value=None, ic=None):
        self.part_id = part_id
        self.value = value
        self.n1 = n1
        self.n2 = n2
        self.ic = ic

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
        print self.get_op_info(ports_v)


class Inductor(Component):
    is_nonlinear = False
    is_symbolic = True

    def __init__(self, part_id='L', n1=None, n2=None, value=None, ic=None):
        self.value = value
        self.n1 = n1
        self.n2 = n2
        self.part_id = part_id
        self.ic = ic
        self.coupling_devices = []


class InductorCoupling(Component):
    is_nonlinear = False
    is_symbolic = True

    def __init__(self, part_id='K', L1=None, L2=None, value_L1=None, value_L2=None):
        self.value_L1 = value_L1
        self.value_21 = value_L2
        self.L1 = L1
        self.L2 = L2

    def __str__(self):
        return "%s %s %g" % (self.L1, self.L2, self.value)

    def get_other_inductor(self, Lselected):
        Lret = None
        if Lselected.upper() == self.L1.upper():
            Lret = self.L2
        elif Lselected.upper() == self.L2.upper():
            Lret = self.L1
        if Lret is None:
            raise Exception, "Mutual inductors bug."
        return Lret

#
# SOURCES
#


class ISource(Component):

    """Generic (ideal) current source:
    Defaults to a DC current source. To implement a time-varying source:
    set _time_function to an appropriate function(time) and is_timedependent=True

    n1: + node
    n2: - node
    dc_value: DC current (A)

    Note: if DC voltage is set and is_timedependent == True, dc_value will be returned
    if the current is evaluated in a DC analysis. This may be useful to simulate a OP
    and then perform a transient analysis with the OP as starting point.
    Otherwise the value in t=0 is used for DC analysis.
    """
    is_nonlinear = False
    is_symbolic = True
    is_timedependent = False
    _time_function = None

    def __init__(self, part_id='I', n1=None, n2=None, dc_value=None, ac_value=0):
        self.part_id = part_id
        self.dc_value = dc_value
        self.abs_ac = numpy.abs(ac_value) if ac_value else None
        self.arg_ac = numpy.angle(ac_value) if ac_value else None
        self.n1 = n1
        self.n2 = n2

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
        """Returns the current in A at the time supplied.
        If time is not supplied, or set to None, or the source is DC, returns dc_value

        This simulator uses Normal convention:
        A positive currents flows in a element from the + node to the - node
        """
        if not self.is_timedependent or (self._time_function == None) or (time == None and self.dc_value is not None):
            return self.dc_value
        else:
            return self._time_function.value(time)


class VSource(Component):

    """Generic (ideal) voltage source:
    Defaults to a DC voltage source. To implement a time-varying source:
    set _time_function to an appropriate function(time) and is_timedependent=True

    n1: + node
    n2: - node
    dc_value: DC voltage (V)

    Note: if DC voltage is set and is_timedependent == True, dc_value will be returned
    if the voltage is evaluated in a DC analysis. This may be useful to simulate a OP
    and then perform a transient analysis with the OP as starting point.
    Otherwise the value in t=0 is used for DC analysis.
    """
    is_nonlinear = False
    is_symbolic = True
    is_timedependent = False
    _time_function = None
    dc_guess = None  # defined in init

    def __init__(self, part_id='V', n1=None, n2=None, dc_value=1.0, ac_value=0):
        self.part_id = part_id
        self.dc_value = dc_value
        self.n1 = n1
        self.n2 = n2
        self.abs_ac = numpy.abs(ac_value) if ac_value else None
        self.arg_ac = numpy.angle(ac_value) if ac_value else None
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
        """Returns the voltage in V at the time supplied.
        If time is not supplied, or set to None, or the source is DC, returns dc_value"""
        if not self.is_timedependent or \
            (self._time_function is None) or \
                (time is None and self.dc_value is not None):
            return self.dc_value
        else:
            return self._time_function.value(time)


class EVSource(Component):

    """Linear voltage controlled voltage source (ideal)

    Source port is a open circuit, dest. port is a ideal voltage source:
    (Vn1 - Vn2) = alpha * (Vsn1 - Vsn2)

    n1: + node, output port
    n2: - node, output port
    sn1: + node, source port
    sn2: - node, source port
    alpha: prop constant between voltages

    """
    is_nonlinear = False
    is_symbolic = True

    def __init__(self, part_id='E', n1=None, n2=None, value=None, sn1=None, sn2=None):
        self.part_id = part_id
        self.n1 = n1
        self.n2 = n2
        self.alpha = value
        self.sn1 = sn1
        self.sn2 = sn2

    def __str__(self):
        return "alpha=%s" % self.alpha


class GISource(Component):

    """Linear voltage controlled current source

    Source port is a short circuit, dest. port is a ideal current source:
    Io = alpha * Is

    Where a positive I enters in n+ and exits from n-

    n1: + node, output port
    n2: - node, output port
    sn1: + node, source port
    sn2: - node, source port
    alpha: prop constant between currents

    """
    is_nonlinear = False
    is_symbolic = True

    def __init__(self, part_id='G', n1=None, n2=None, value=None, sn1=None, sn2=None):
        self.part_id = part_id
        self.n1 = n1
        self.n2 = n2
        self.alpha = value
        self.sn1 = sn1
        self.sn2 = sn2

    def __str__(self):
        return "value=%s" % self.alpha


class HVSource(Component):  # TODO: fixme

    """Linear current controlled voltage source
    """
    is_nonlinear = False
    is_symbolic = True

    def __init__(self, part_id='H', n1=None, n2=None, value=None, sn1=None, sn2=None):
        print "HVSource not implemented. TODO"
        self.part_id = part_id
        self.n1 = n1
        self.n2 = n2
        self.alpha = value
        self.sn1 = sn1
        self.sn2 = sn2
        import sys
        sys.exit(1)

    def __str__(self):
        raise Exception, "HVSource not implemented. TODO"


#
# Functions for time dependent sources  #
#

class pulse:
    # PULSE(V1 V2 TD TR TF PW PER)
    _type = "V"
    v1 = None
    v2 = None
    td = None
    per = None
    tr = None
    tf = None
    pw = None

    def __init__(self, v1=None, v2=None, td=None, tr=None, pw=None, tf=None, per=None):
        self.v1 = v1
        self.v2 = v2
        self.td = td
        self.per = per
        self.tr = tr
        self.tf = tf
        self.pw = pw

    def value(self, time):
        if not self.ready():
            print "Error: pulse function not well defined. This is a bug."
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

    def ready(self):
        if self.v1 == None or self.v2 == None or self.td == None or self.tr == None or self.pw == None or \
                self.tf == None or self.per == None:
            return False
        else:
            return True

    def __str__(self):
        return "type=pulse " + \
            self._type.lower() + "1=" + str(self.v1) + " " + \
            self._type.lower() + "2=" + str(self.v2) + \
            " td=" + str(self.td) + " per=" + str(self.per) + \
            " tr=" + str(self.tr) + " tf=" + str(self.tf) + \
            " pw=" + str(self.pw)


class sin:
    # SIN(VO VA FREQ TD THETA)
    td = None
    vo = None
    va = None
    freq = None
    theta = None
    _type = "V"

    def __init__(self, vo=None, va=None, freq=None, td=None, theta=None):
        self.vo = vo
        self.va = va
        self.freq = freq
        self.td = td
        self.theta = theta

    def value(self, time):
        if not self.ready():
            printing.print_general_error(
                "Error: sin function not well defined. This is a bug.")
        if time < self.td:
            return self.vo
        elif self.theta:
            return self.vo + self.va * math.exp(-1 * (time - self.td) / self.theta) * math.sin(2 * math.pi * self.freq * (time - self.td))
        else:
            return self.vo + self.va * math.sin(2 * math.pi * self.freq * (time - self.td))

    def ready(self):
        if self.vo == None or self.va == None or self.freq == None or self.td == None or self.theta == None:
            return False
        else:
            return True

    def __str__(self):
        return "type=sin " + \
            self._type.lower() + "o=" + str(self.vo) + " " + \
            self._type.lower() + "a=" + str(self.va) + \
            " freq=" + str(self.freq) + " theta=" + str(self.theta) + \
            " td=" + str(self.td)


class exp:
    # EXP(V1 V2 TD1 TAU1 TD2 TAU2)
    v1 = None
    v2 = None
    td1 = None
    tau1 = None
    td2 = None
    tau2 = None
    _type = "V"

    def __init__(self, v1=None, v2=None, td1=None, tau1=None, td2=None, tau2=None):
        self.v1 = v1
        self.v2 = v2
        self.td1 = td1
        self.tau1 = tau1
        self.td2 = td2
        self.tau2 = tau2

    def value(self, time):
        if not self.ready():
            printing.print_general_error(
                "Error: exp function not well defined. This is a bug.")
        if time < self.td1:
            return self.v1
        elif time < self.td2:
            return self.v1 + (self.v2 - self.v1) * (1 - math.exp(-1 * (time - self.td1) / self.tau1))
        else:
            return self.v1 + (self.v2 - self.v1) * (1 - math.exp(-1 * (time - self.td1) / self.tau1)) + (self.v1 - self.v2) * (1 - math.exp(-1 * (time - self.td2) / self.tau2))

    def ready(self):
        if self.v1 == None or self.v2 == None or self.td1 == None or self.tau1 == None or self.td2 == None \
                or self.tau2 == None:
            return False
        return True

    def __str__(self):
        return "type=exp " + \
            self._type.lower() + "1=" + str(self.v1) + " " + \
            self._type.lower() + "2=" + str(self.v2) + \
            " td1=" + str(self.td1) + " td2=" + str(self.td2) + \
            " tau1=" + str(self.tau1) + " tau2=" + str(self.tau2)
