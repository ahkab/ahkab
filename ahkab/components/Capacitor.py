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
from .Component import Component
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
