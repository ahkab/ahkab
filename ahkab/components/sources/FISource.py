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

