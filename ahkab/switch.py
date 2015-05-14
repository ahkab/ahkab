# -*- coding: iso-8859-1 -*-
# switch.py
# Implementation of the voltage controlled switch model
# Copyright 2013 Giuseppe Venturini
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
Implementation of a voltage controlled switch.

This module defines two classes: switch_device, switch_model

"""

# sn1 o--+         +--o n1
#        |         |
#       +-+      \ o
#       |R|       \
#       +-+        +
#        |         |
# sn2 o--+         +--o n2

from __future__ import (unicode_literals, absolute_import,
                        division, print_function)

import math

from . import options
from . import printing


class switch_device:

    """This is a general switch element.

    It has the following structure:

    .. image:: images/elem/switch1.svg

    |

    In ASCII for those who are consulting the documentation from the
    Python command line:

    ::

        sn1 o--+         +--o n1
               |         |
              +-+      \ o
              |R|       \\
              +-+        +
               |         |
        sn2 o--+         +--o n2

    The behavior is set by the model supplied.

    The device instance calls the following methods in the model:

    * ``get_i(ports_v, device)`` - output current
    * ``get_go(ports_v, device)`` - ouput conductance
    * ``get_gm(ports_v, device)`` - output transconductance
    * ``get_dc_guess(self, is_on)`` - guesses for OP

    The device instance accesses the following attributes:
    ``part_id`` (a string), the device label.

    """

    def __init__(self, n1, n2, sn1, sn2, model, ic=None, part_id='S'):
        """
        **Parameters:**

        n1 : str
            Positive output node (+)
        n2 : str
            Negative output node (-)
        sn1 : str
            Positive input node (+)
        sn2 : str
            Negative input node (-)
        model : model obj
            An instance of (v)switch_model
        ic : bool, optional
            The initial conditions: ``True`` stands for on, ``False`` for off.

        Selected methods:

        - :func:`get_output_ports` -> (n1, n2)
        - :func:`get_drive_ports` -> (n1, n2), (ns1, ns2)

        """
        class dev_class:
            pass
        self.device = dev_class()
        self.device.is_on = ic if ic is not None else False
        self.sn1 = sn1
        self.sn2 = sn2
        self.n1 = n1
        self.n2 = n2
        self.ports = ((self.n1, self.n2), (self.sn1, self.sn2))
        self.model = model
        self.opdict = {}
        self.opdict.update({'state': (float('nan'), float('nan'))})
        self.part_id = part_id
        self.is_nonlinear = True
        self.is_symbolic = True
        self.dc_guess = self.model.get_dc_guess(self.device.is_on)

    def get_drive_ports(self, op):
        """Get the ports that drive the output ports.

        **Parameters:**

        op : op solution
            The OP where the drive ports are used.

        **Returns:**

        pts : tuple of tuples of ports nodes, as: ``(port0, port1, port2 ... )``

        Where each port is in the form: ``port0 = (nplus, nminus)``
        """
        return self.ports

    def get_output_ports(self):
        """Get the output port.

        The output port is ``(n1, n2)`` for the voltage-controlled switch case.

        **Returns:**

        pts : tuple of tuples of ports nodes
            Such as: ``(port0, port1, port2 ... )``.
            Where each port is in the form: ``port0 = (nplus, nminus)``
        """
        return ((self.n1, self.n2),)

    def __str__(self):
        rep = self.model.name + " " + str(self.device.is_on)
        return rep

    def i(self, op_index, ports_v, time=0):
        """Returns the current flowing in the element.

        The element is assumed to be biased with the voltages
        applied as specified in the ``ports_v`` vector.

        **Parameters:**

        op_index : int
            The index of the output port for which the current is evaluated.
        ports_v : tuple
            A tuple constructed such as ``(voltage_across_port0, voltage_across_port1, ... )``
        time : float, optional
            The simulation time at which the evaluation is performed. It is
            needed by time-variant elements, and it has no effect here. Set it
            to ``None`` during DC analysis.

        **Returns:**

        i : int
            The output current.
        """
        ret = self.model.get_i(ports_v, self.device)
        # This may be used for debugging
        # print str(ports_v)+" Isw: %g\tRo: %g\tgm: %g" % (ret, 1/self.g(0,
        # ports_v, 0), self.g(0, ports_v, 1))
        return ret

    def update_status_dictionary(self, ports_v):
        """Updates an internal dictionary that can then be used to provide
        information to the user regarding the status of the element.

        Normally, one would call :func:`get_op_info`.

        **Returns:**

        ``None``.
        """
        if self.opdict is None:
            self.opdict = {}
        if not (self.opdict['state'] == ports_v[0] and 'R' in self.opdict):
            self.opdict['state'] = ports_v[0]
            self.opdict['R'] = float(1.0 / self.g(0, ports_v[0], 0))
            self.opdict['I'] = float(self.i(0, ports_v[0]))
            self.opdict['STATUS'] = self.device.is_on

    def get_op_info(self, ports_v):
        """Information regarding the Operating Point (OP)

        **Parameters:**

        ports_v : list of lists
            The parameter is to be set to ``[[v]]``, where ``v`` is the voltage
            applied to the switch terminals.

        **Returns:**

        op_keys : list of strings
            The labels corresponding to the numeric values in ``op_info``.
        op_info : list of floats
            The values corresponding to ``op_keys``.
        """
        self.update_status_dictionary(ports_v)
        status = "ON" if self.opdict['STATUS'] else "OFF"
        op_keys = ['Part ID', 'STATUS', "VO [V]", "VS [V]", u"R [\u2126]",
                   "I [A]"]
        op_info = [self.part_id, status, float(self.opdict['state'][0]),
                   float(self.opdict['state'][1]), self.opdict["R"],
                   self.opdict['I']]
        return op_keys, op_info

    def g(self, op_index, ports_v, port_index, time=0):
        """Returns the differential (trans)conductance.

        The transconductance is computed wrt the port specified by
        ``port_index`` when the element has the voltages specified in
        ``ports_v`` across its ports, at (simulation) ``time``.

        **Parameters:**

        ports_v : list
            Voltages applied to the switch. The list should be in the form:
            ``[voltage_across_port0, voltage_across_port1, ... ]``
        port_index : int
            The index of the output port.
        time : float
            The simulation time at which the evaluation is performed. Set it to
            ``None`` during DC analysis.

        **Returns:**

        g : float
            The transconductance.
        """

        assert op_index == 0
        assert port_index < 2

        if port_index == 0:
            return self.model.get_go(ports_v, self.device)
        if port_index == 1:
            return self.model.get_gm(ports_v, self.device)
        else:
            raise Exception("Unknown port index passed to switch: bug")

    def get_value_function(self, identifier):
        def get_value(self):
            return self.opdict[identifier]
        return get_value

    def get_netlist_elem_line(self, nodes_dict):
        """Return a netlist line corresponding to the switch."""
        return "%s %s %s %s %s %s %s" % (self.part_id, nodes_dict[self.n1],
                                nodes_dict[self.n2], nodes_dict[self.sn1],
                                nodes_dict[self.sn2], self.model.name, \
                                str(self.device.is_on))


VT_DEFAULT = 0.0
VH_DEFAULT = 0.0
RON_DEFAULT = 1.
ROFF_DEFAULT = 1. / options.gmin


class vswitch_model:
    """Voltage-controlled switch model.

    ::

        sn1 o--+         +--o n1
               |         |
              +-+      \ o
              |R|       \\
              +-+        +
               |         |
        sn2 o--+         +--o n2


    Note that:

    * R is infinite.
    * The voltage needed to close the switch is:
      :math:`V(s_{n1})-V(s_{n2}) > V_T+V_H`.
    * To re-open it, one needs to satisfy the relationship:
      :math:`V(s_{n1})-V(s_{n2}) < V_T-V_H`.

    The switch commutes between two statuses:

    * :math:`R_{OUT} = R_{OFF}`
    * :math:`R_{OUT} = R_{ON}`

    None of which can be set to zero or infinite.

    The switching characteristics are modeled with :math:`tanh(x)`.
    """

    def __init__(self, name, VT=None, VH=None, VON=None, VOFF=None, RON=None, ROFF=None):
        self.name = name
        # convert to VT and VH
        if VON is not None or VOFF is not None:
            VT, VH = self._get_VTVH_from_VONVOFF(float(VON), float(VOFF))
        self.VT = float(VT) if VT is not None else VT_DEFAULT
        self.VH = float(VH) if VH is not None else VH_DEFAULT
        self.RON = float(RON) if RON is not None else RON_DEFAULT
        self.ROFF = float(ROFF) if ROFF is not None else ROFF_DEFAULT
        self.A = (self.RON - self.ROFF) / 2
        self.B = (self.RON + self.ROFF) / 2.
        self.is_on = False
        self._set_status(self.is_on)
        self.SLOPE = 1e2

    def _get_VTVH_from_VONVOFF(self, VON, VOFF):
        if VON is None or VOFF is None:
            raise ValueError
        VT = (VON - VOFF) / 2.0 + VOFF
        return VT, VT * 1e-3

    def _get_V(self, is_on):
        """Get the effective switching voltage (hyst taken into account)
        """
        return self.VT + self.VH * 2 * (not is_on) - self.VH

    def _set_status(self, is_on):
        """Set the switch status, which meeans setting the effective
        switching voltage self.V (w hyst taken into account)
        """
        self.V = self.VT + self.VH * 2 * (not is_on) - self.VH

    def _update_status(self, vin, dev, debug=False):
        """Check the switch status and move to the other if needed.
        """
        Vtest = self._get_V(dev.is_on)
        R1 = self.A * math.tanh((vin - Vtest) * self.SLOPE) + self.B
        Vtest = self._get_V(not dev.is_on)
        R2 = self.A * math.tanh((vin - Vtest) * self.SLOPE) + self.B
        self._set_status(dev.is_on)
        if vin > self.V and not dev.is_on and R1 - R2 == 0.0:
            if debug:
                print("Switching ON: %g" % (vin,))
            dev.is_on = True
            self._set_status(dev.is_on)
        if vin < self.V and dev.is_on and R1 - R2 == 0.0:
            if debug:
                print("Switching OFF: %g" % (vin,))
            dev.is_on = False
            self._set_status(dev.is_on)
        self.is_on = dev.is_on

    def get_dc_guess(self, is_on):
        """Returns a list of two floats to be used as initial guesses for the OP analysis
        """
        return [self.VT * (.9 + is_on * .2)] * 2

    def print_model(self):
        """All the internal parameters of the model get printed out,
        for visual inspection.
        """
        arr = []
        arr.append(
            [self.name, "", "", "SWITCH MODEL", "", "", "", "", "",  "", "", ""])
        arr.append(["VT", "[V]", self.VT, "VH", "[V]:", self.VH,
                   "RON", "[ohm]", self.RON, "ROFF", "[ohm]", self.ROFF])
        print(printing.table(arr))

    def get_i(self, xxx_todo_changeme, dev, debug=False):
        """Returns the output current.
        """
        (vout, vin) = xxx_todo_changeme
        self._update_status(vin, dev)
        R = self.A * math.tanh((vin - self.V) * self.SLOPE) + self.B
        return vout / R

    def get_go(self, xxx_todo_changeme1, dev, debug=False):
        """Returns the output conductance d(I)/d(Vn1-Vn2)."""
        (vout, vin) = xxx_todo_changeme1
        self._update_status(vin, dev)
        R = self.A * math.tanh((vin - self.V) * self.SLOPE) + self.B
        return 1. / R

    def get_gm(self, xxx_todo_changeme2, dev, debug=False):
        """Returns the source to output transconductance or d(I)/d(Vsn1-Vsn2)."""
        (vout, vin) = xxx_todo_changeme2
        self._update_status(vin, dev)
        gm = self.A * self.SLOPE * (math.tanh(self.SLOPE * (self.V - vin)) ** 2 - 1) / (
            self.A * math.tanh(self.SLOPE * (self.V - vin)) - self.B) ** 2
        return gm + options.gmin

