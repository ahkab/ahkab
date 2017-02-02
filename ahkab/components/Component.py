# -*- coding: iso-8859-1 -*-
# Component.py
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

Defining new elements and subclassing ``Component``
---------------------------------------------------

We recommend to subclass :class:`ahkab.components.Component` if you intend to
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
