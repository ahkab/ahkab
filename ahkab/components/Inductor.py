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

