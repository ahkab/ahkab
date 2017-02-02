# -*- coding: iso-8859-1 -*-
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
from ..Component import Component
import numpy as np

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

