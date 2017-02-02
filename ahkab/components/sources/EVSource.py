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

