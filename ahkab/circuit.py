# -*- coding: iso-8859-1 -*-
# circuit.py
# Describes the circuit
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
Introduction
\"\"\"\"\"\"\"\"\"\"\"\"

A circuit is described in the ahkab simulator by an instance of the
:class:`Circuit` class.

This class holds everything is needed to simulate the circuit, except
the specification of the analyses to be performed.

To rewrite a netlist from a Circuit instance see the
:mod:`ahkab.printing` module.

The Circuit
\"\"\"\"\"\"\"\"\"\"\"

A circuit is derived from a list which contains all its elements.

Conceptually, every time an element is to be inserted in the circuit,
two operations have to be performed:

* The element must be appended to the ``Circuit`` instance.
* Its connections should be ensure checking that
  the nodes the element refers to are indeed existing circuit nodes.

To simplify the operation of adding a component to a ``Circuit``,
the following convenience methods are provided to the user to add and
remove most elements to the circuit:

* :func:`Circuit.add_resistor`
* :func:`Circuit.add_capacitor`
* :func:`Circuit.add_inductor`
* :func:`Circuit.add_vsource`
* :func:`Circuit.add_isource`
* :func:`Circuit.add_diode`
* :func:`Circuit.add_mos`
* :func:`Circuit.add_cccs`
* :func:`Circuit.add_vcvs`
* :func:`Circuit.add_vccs`
* :func:`Circuit.add_user_defined`
* :func:`Circuit.remove_elem`

Example:

.. code-block:: python

    mycircuit = circuit.Circuit(title="Example circuit", filename=None)
    # no filename since there will be no deck associated with this circuit.
    # get the ref node (gnd)
    gnd = mycircuit.get_ground_node()
    # add a node named n1 and a 600 ohm resistor connected between n1 and gnd
    mycircuit.add_resistor(part_id="R1", n1="n1", n2=gnd, R=600)

Refer to the methods help for additional information.

Nodes
\"\"\"\"\"

The nodes are internally stored in the following way: we assign to each node an
internal ID, independetly from its external identifier used in the netlist.
Those IDs are integers.

The simulator uses always the internal names. When the results are
presented to the user, the internal node is not showed, the external
identifier (or external node name) is printed instead.

This is done through:

.. code-block:: python

    my_circuit = Circuit()
    ...
    [ init code ]
    ...
    print "This is a node" + my_circuit.nodes_dict[int_node]

.. rubric:: Internal only nodes

The number of internal only nodes (added automatically by the simulator)
is held in ``Circuit.internal_nodes``. That value shouldn't be changed by
hand.

Device models
\"\"\"\"\"\"\"\"\"\"\"\"\"

Non-linear elements have their operation described by specialized routines
held in their module.

They are stored in ``Circuit.models`` (of type dict), the following methods
are provided to add and remove device models to a Circuit instance.

* :func:`Circuit.add_model`
* :func:`Circuit.remove_model`

Reference
\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\""\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"

"""

from __future__ import (unicode_literals, absolute_import,
                        division, print_function)

import importlib
import importlib.util

import copy
import math

from . import components
from . import diode
from . import ekv
from . import mosq
from . import printing
from . import py3compat
from . import switch

# will be added here by netlist_parser and circuit instances
user_defined_modules_dict = {}

class Circuit(list):
    """The circuit class.

    **Parameters:**

    title : string
        The circuit title.

    filename : string, optional

        .. deprecated:: 0.09

        If the circuit instance corresponds to a netlist file on disk,
        set this to the netlist filename.

    """
    def __init__(self, title, filename=None):
        self.title = title
        self.filename = filename
        self.nodes_dict = {}  # {int_node:ext_node, int_node:ext_node}
        self.internal_nodes = 0
        self.models = {}
        self.gnd = '0'

    def __str__(self):
        s = "* " + self.title + "\n"
        for elem in self:
            s += elem.get_netlist_elem_line(self.nodes_dict) + "\n"
        return s[:-1]

    def create_node(self, name):
        """Creates a new circuit node

        If there is a node in the circuit with the same name, ValueError is
        raised.

        **Parameters:**

        name : string
            the _unique_ identifier of the node.

        **Returns:**

        node : string
            the _unique_ identifier of the node, to be used for subsequent
            element declarations, for example.

        :raises ValueError: if a new node with the given id cannot be created,
          for example because a node with the same name already exists in the
          circuit. The only exception is the ground node, which has the
          reserved id ``'0'``, and for which this method won't raise any
          exception.
        :raises TypeError: if the parameter ``name`` is not of "text" type (what
          that means exactly depends on which version of Python you are using.)

        """
        if type(name) not in py3compat.string_types:
            raise TypeError("The node %s should have been of text type" %
                            name)
        got_ref = 0 in self.nodes_dict
        if name not in self.nodes_dict:
            if name == '0':
                int_node = 0
            else:
                int_node = int(len(self.nodes_dict)/2) + 1*(not got_ref)
            self.nodes_dict.update({int_node:name})
            self.nodes_dict.update({name:int_node})
        else:
            raise ValueError('Impossible to create new node %s: node exists!'
                             % name)
        return name

    def add_node(self, ext_name):
        """Adds the supplied node to the circuit, if needed.

        When a 'normal' (not the reference) node is added, a internal
        name (or label) is assigned to it.

        The nodes labels are stored in ``Circuit.nodes_dict``, as a dictionary of pairs
        like ``{int_node:ext_node}``.

        Those internal names are integers, by definition, and are generated
        starting from 1, then 2, 3, 4, 5...
        The integer ``0`` is reserved for the reference node (gnd), which is required
        for the circuit to be non-pathological and has ``ext_name=str(int_name)='0'``.

        Notice that this method doesn't halt or print errors if the node is already been
        added previsiously. It simply returns the internal node name assigned to it.

        **Parameters:**

        ext_name : string
            The unique identifier of the node.

        **Returns:**

        int_name : string
            the *unique* *internal* ciecuit identifier of the node.

        :raises TypeError: if the parameter ``ext_name`` is not of "text" type
          (what that means exactly depends on which version of Python you are
          using.)
        """
        # must be text (str unicode...)
        if type(ext_name) not in py3compat.string_types:
            raise TypeError("The node %s should have been of text type" %
                            ext_name)
        # test: do we already have it in the dictionary?
        if ext_name not in self.nodes_dict:
            if ext_name == '0':
                int_node = 0
            else:
                got_ref = 0 in self.nodes_dict
                int_node = int(len(self.nodes_dict)/2) + 1*(not got_ref)
            self.nodes_dict.update({int_node:ext_name})
            self.nodes_dict.update({ext_name:int_node})
        else:
            int_node = self.nodes_dict[ext_name]
        return int_node

    def new_internal_node(self):
        """Generate implicit internal nodes.

        Some devices are made of a group of other devices, connected by
        "internal only" nodes, which have the prefix ``'INT'`` and the
        simulator treats specially, hiding them from the user if not
        explicitly asked about them.

        This method generates the external names for such nodes and inserts them
        in the circuit.

        **Returns:**

        ext_node : string
            The corresponding external node id.
        """

        ext_node = "INT" + str(self.internal_nodes)
        self.internal_nodes = self.internal_nodes + 1
        self.create_node(ext_node)
        return ext_node

    def get_nodes_number(self):
        """Returns the number of nodes in the circuit"""
        return int(len(self.nodes_dict)/2)

    def is_int_node_internal_only(self, int_node):
        """Check whether an internal node is an "internal only node" or not.

        **Parameters:**

        int_node : int
            The internal only node to be checked.

        **Returns:**

        chk : boolean
            The result of the check.

        :raises TypeError: if the supplied node is not an ``int``. Typically
          this happens when the method is called with an *external* name.
        """
        if type(int_node) is not int:
            raise TypeError('Expecting an INTERNAL node of type int, got %s.' %
                            type(int_node))
        return self.nodes_dict[int_node].find("INT") > -1

    def is_nonlinear(self):
        """Check whether the circuit is non-linear or not.

        **Returns:**

        chk : boolean
            The result of the check.
        """
        for elem in self:
            if elem.is_nonlinear:
                return True
        return False

    def get_locked_nodes(self):
        """Get all nodes connected to non-linear elements.

        This list is meant to be passed to ``dc_solve`` or ``mdn_solver`` to be
        used in ``get_td`` to evaluate the damping coefficient in a
        Newton-Rhapson iteration.

        **Returns:**

        locked_nodes : list
            A list of internal nodes.
        """
        locked_nodes = []
        nl_elements = [elem for elem in self if elem.is_nonlinear]
        for elem in nl_elements:
            oports = elem.get_output_ports()
            for index in range(len(oports)):
                ports = elem.get_drive_ports(index)
                for port in ports:
                    locked_nodes.append(port)
        return locked_nodes

    def ext_node_to_int(self, ext_node):
        """This function returns the integer id associated with an external node id.

        **Parameters:**

        ext_node : string
            The external node id to be converted.

        **Returns:**

        int_node : int
            The internal node associated.

        """
        return self.nodes_dict[ext_node]

    def int_node_to_ext(self, int_node):
        """This function returns the string id associated with the integer internal node id
        ``int_node``.

        **Parameters:**

        int_node : int
            The internal node id to be converted.

        **Returns:**

        ext_node : string
            the string id associated with ``int_node``.
        """
        return self.nodes_dict[int_node]

    def has_duplicate_elem(self):
        """Self-check for duplicate elements.

        No circuit should ever have duplicate elements
        (ie elements with the same ``part_id``).

        **Returns:**

        chk : boolean
            The result of the check.
        """
        all_ids = tuple(map(lambda e: e.part_id, self))
        return len(set(all_ids)) != len(all_ids)

    def get_ground_node(self):
        """Returns the reference node, AKA GND."""
        return '0'

    def get_elem_by_name(self, part_id):
        """Get a circuit element from its ``part_id`` value.

        If no matching element is found, the method returns
        ``None``. This may change in the future.

        **Parameters:**

        part_id : string
            The ``part_id`` of the element

        **Returns:**

        elem : circuit element
            Depending whether a matching element was found or not.

        :raises ValueError: if the element is not found.
        """
        for e in self:
            if e.part_id.lower() == part_id.lower():
                return e
        raise ValueError('Element %s not found' % part_id)

    def add_model(self, model_type, model_label, model_parameters):
        """Add a model to the available circuit models.

        **Parameters:**

        model_type : string
            the model type (eg "ekv"). Right now, the possible values are:
            ``"mosq"``, ``"ekv"``, ``"diode"``, ``"sw"``.

        model_label : string
            a unique identifier for the model being added (eg. ``"nch1"``).

        model_parameters: dict
            a dictionary holding the parameters to be supplied to the
            model to instantiate it.

        """

        if 'name' not in model_parameters:
            model_parameters.update({'name':model_label})
        if model_type == "ekv":
            model_iter = ekv.ekv_mos_model(**model_parameters)
            model_iter.name = model_label
        elif model_type == "mosq":
            model_iter = mosq.mosq_mos_model(**model_parameters)
            model_iter.name = model_label
        elif model_type == "diode":
            model_iter = diode.diode_model(**model_parameters)
            model_iter.name = model_label
        elif model_type == "sw":
            model_iter = switch.vswitch_model(**model_parameters)
            model_iter.name = model_label
        else:
            raise CircuitError("Unknown model type %s" % (model_type,))
        self.models.update({model_label: model_iter})

    def remove_model(self, model_label):
        """Remove a model from the available models.

        **Parameters:**

        model_label : string
            the unique identifier corresponding to the model
            being removed.

        .. note::

            This method currently silently ignores models that are not defined.

        """
        if self.models is not None and model_label in self.models:
            del self.models[model_label]
        # should print a warning here

    def add_resistor(self, part_id, n1, n2, value):
        """Adds a resistor to the circuit.

        The resistor instance is added to the circuit elements
        and connected to the provided nodes. If the nodes are not
        found in the circuit, they are created and added as well.

        **Parameters:**

        part_id : string
            the resistor part_id (eg "R1"). The first letter is replaced by an R

        n1, n2 : string
            the nodes to which the resistor is connected.

        value : float,
            The resistance between ``n1`` and ``n2`` in Ohm.

        .. seealso::

            :func:`add_resistor`, :func:`add_capacitor`,
            :func:`add_inductor`, :func:`add_vsource`, :func:`add_isource`,
            :func:`add_diode`, :func:`add_mos`, :func:`add_vcvs`, :func:`add_vccs`,
            :func:`add_cccs`, :func:`add_user_defined`, :func:`remove_elem`

        """
        n1 = self.add_node(n1)
        n2 = self.add_node(n2)

        if value == 0:
            raise CircuitError("ZERO-valued resistors are not allowed.")

        elem = components.Resistor(part_id=part_id, n1=n1, n2=n2, value=value)
        self.append(elem)

    def add_capacitor(self, part_id, n1, n2, value, ic=None):
        """Adds a capacitor to the circuit.

        The capacitor instance is added to the circuit elements
        and connected to the provided nodes. If the nodes are not
        found in the circuit, they are created and added as well.

        **Parameters:**

        part_id : string
            The capacitor part_id (eg "C1"). The first letter is always C.

        n1, n2 : string
            The nodes to which the element is connected.

        value : float
            The capacitance value.

        ic : float, optional
            The initial condition, if any. See the simulation docs for
            how this affects the results.

        .. seealso::
            :func:`add_resistor`,
            :func:`add_inductor`, :func:`add_vsource`, :func:`add_isource`,
            :func:`add_diode`, :func:`add_mos`, :func:`add_vcvs`, :func:`add_vccs`,
            :func:`add_cccs`, :func:`add_user_defined`, :func:`remove_elem`

        """
        if value == 0:
            raise CircuitError("ZERO-valued capacitors are not allowed.")

        n1 = self.add_node(n1)
        n2 = self.add_node(n2)

        elem = components.Capacitor(part_id=part_id, n1=n1, n2=n2, value=value, ic=ic)

        self.append(elem)

    def add_inductor(self, part_id, n1, n2, value, ic=None):
        """Adds an inductor to the circuit.

        The inductor instance is added to the circuit elements
        and connected to the provided nodes. If the nodes are not
        found in the circuit, they are created and added as well.

        **Parameters:**

        part_id : string
            The inductor part_id (eg "Lfilter"). The first letter is always L.

        n1, n2 : string
            The nodes to which the element is connected. Eg. ``"in"`` or ``"out_a"``.

        value : float
            The inductance value.

        ic : float, optional
            Initial condition, see simulation types for how this affects
            the results.

        .. seealso::

            :func:`add_resistor`, :func:`add_capacitor`,
            :func:`add_inductor`, :func:`add_vsource`, :func:`add_isource`,
            :func:`add_diode`, :func:`add_mos`, :func:`add_vcvs`, :func:`add_vccs`,
            :func:`add_cccs`, :func:`add_user_defined`, :func:`remove_elem`
        """

        n1 = self.add_node(n1)
        n2 = self.add_node(n2)

        elem = components.Inductor(part_id=part_id, n1=n1, n2=n2, value=value, ic=ic)

        self.append(elem)

    def add_inductor_coupling(self, part_id, L1, L2, value):
        """Add a coupling between two inductors.

        **Parameters:**

        part_id : string
            The part ID for the inductor coupling device. Eg. ``'K1'``,
            the first letter is always ``'K'``.
        L1 : string
            The part ID of the first inductor to be coupled.
        L2 : string
            The part ID of the second inductor to be coupled.
        value : float
            The ``k`` value of the mutual coupling coefficient.
            Its value must be greater than zero and lesser or equal to``1``
            or instability ensues.
        """
        L1elem, L2elem = None, None

        for e in self:
            if isinstance(e, components.Inductor) and (L1 == e.part_id):
                L1elem = e
            elif isinstance(e, components.Inductor) and (L2 == e.part_id):
                L2elem = e

        if L1elem is None or L2elem is None:
            error_msg = "One or more coupled inductors for %s were not found: %s (found: %s), %s (found: %s)." % \
                (part_id, L1, L1elem is not None, L2, L2elem is not None)
            raise ValueError(error_msg)

        M = math.sqrt(L1elem.value * L2elem.value) * value

        elem = components.InductorCoupling(part_id=part_id, L1=L1, L2=L2, K=value, M=M)
        L1elem.coupling_devices.append(elem)
        L2elem.coupling_devices.append(elem)

        self.append(elem)

    def add_vsource(self, part_id, n1, n2, dc_value, ac_value=0, function=None):
        """Adds a voltage source to the circuit (also takes care that the nodes
        are added as well).

        **Parameters:**

        part_id : string
            The voltage source part_id (eg "VA"). The first letter is always V.
        n1, n2 : string
            The nodes to which the element is connected. Eg. ``"in"`` or
            ``"out_a"``.
        dc_value : float
            DC voltage value
        ac_value : float, optional
            AC voltage value, defaults to 0.
        function : function, optional
            Time function. See devices.py for built-in options.
        """
        n1 = self.add_node(n1)
        n2 = self.add_node(n2)

        elem = components.sources.VSource(part_id=part_id, n1=n1, n2=n2, dc_value=dc_value,
                               ac_value=ac_value)

        if function is not None:
            elem.is_timedependent = True
            elem._time_function = function

        self.append(elem)

    def add_isource(self, part_id, n1, n2, dc_value, ac_value=0, function=None):
        """Adds a current source to the circuit (also takes care that the nodes
        are added as well).

        **Parameters:**

        part_id : string
            The current source ID (eg ``"IA"`` or ``"I3"``). The first letter
            is always I.
        n1, n2 : string
            The nodes to which the element is connected, eg. ``"in"`` or ``"out1"``.
        dc_value : float
            DC current value.
        ac_value :float, optional
            AC current value, defaults to 0.
        function : function, optional
            Time function. See devices.py for built-in options.
        """
        n1 = self.add_node(n1)
        n2 = self.add_node(n2)

        elem = components.sources.ISource(part_id=part_id, n1=n1, n2=n2, dc_value=dc_value,
                               ac_value=ac_value)

        if function is not None:
            elem.is_timedependent = True
            elem._time_function = function

        self.append(elem)

    def add_diode(self, part_id, n1, n2, model_label, models=None, Area=None,
                  T=None, ic=None, off=False):
        """Adds a diode to the circuit (also takes care that the nodes
        are added as well).

        **Parameters:**

        part_id : string
            The diode ID (eg "D1"). The first letter is always D.
        n1, n2 : string
            the nodes to which the element is connected. eg. ``"in"`` or
            ``"out_a"``
        model_label : string
            The diode model identifier. The model needs to be added
            first, then the elements using it.
        models : dict, optional
            List of available model instances. If not set or ``None``,
            the circuit models will be used (recommended).
        Area : float, optional
            Scaled device area (optional, defaults to 1.0)
        T : float, optional
            Operating temperature (no temperature dependence yet)
        ic : float, optional
            Initial condition (not really implemented yet)
        off : bool, optional
            Consider the diode to be initially off.
        """
        n1 = self.add_node(n1)
        n2 = self.add_node(n2)
        if models is None:
            models = self.models
        if model_label not in models:
            raise ModelError("Unknown diode model id: " + model_label)

        elem = diode.diode(part_id=part_id, n1=n1, n2=n2, model=models[
                           model_label], AREA=Area, T=T, ic=ic, off=off)
        self.append(elem)


    def add_mos(self, part_id, nd, ng, ns, nb, w, l, model_label, models=None,
                m=1, n=1):
        """Adds a mosfet to the circuit (also takes care that the nodes
        are added as well).

        **Parameters:**

        part_id : string
            The mos part_id (eg "M1"). The first letter is always M.
        nd : string
            The drain node.
        ng : string
            The gate node.
        ns : string
            The source node.
        nb : string
            The bulk node.
        w : float
            The gate width.
        l : float
            The gate length.
        model_label : string
            The model identifier.
        models : dict, optional
            The circuit models.
        m : int, optional
            Shunt multiplier value. Defaults to 1.
        n : int, optional
            Series multiplier value, not always supported. Defaults to 1.
        """
        nd = self.add_node(nd)
        ng = self.add_node(ng)
        ns = self.add_node(ns)
        nb = self.add_node(nb)

        if models is None:
            models = self.models

        if model_label not in models:
            raise ModelError("Unknown model id: " + model_label)

        if isinstance(models[model_label], ekv.ekv_mos_model):
            elem = ekv.ekv_device(part_id, nd, ng, ns, nb, w, l,
                                  models[model_label], m, n)
        elif isinstance(models[model_label], mosq.mosq_mos_model):
            elem = mosq.mosq_device(part_id, nd, ng, ns, nb, w, l,
                                    models[model_label], m, n)
        else:
            raise Exception("Unknown model type for " + model_label)

        self.append(elem)


    def add_cccs(self, part_id, n1, n2, source_id, value):
        """Adds a current-controlled current source (CCCS) to the circuit

        This method takes care that its nodes are added as well.

        **Parameters:**

        part_id : string
            The cccs ID (eg ``'F1'``). The first letter is always ``'F'``.
        n1, n2 : strings
            The output port nodes, where the output current is
            forced. Eg. "outp", "outm" or "out_a", "out_b".
        source_id : string
            The voltage source to be used to sense the current that drives
            the output. Eg. ``'V1'``.
        value : float
            The proportionality factor between input (:math:`I_s`) and output
            (:math:`I_o`) currents. Mathematically:

            .. math::

                I_o = \\alpha I_s

        .. seealso::

            :class:`ahkab.components.sources.FISource`

        """
        # Add the nodes, this is SAFE: if a node is already known to the circuit,
        # the methods will just ignore the request.
        n1 = self.add_node(n1)
        n2 = self.add_node(n2)
        # instantiate the element
        elem = components.sources.FISource(part_id=part_id, n1=n1, n2=n2,
                                source_id=source_id, value=value)
        # add it!
        self.append(elem)

    def add_ccvs(self, part_id, n1, n2, source_id, value):
        """Adds a current-controlled voltage source (CCCS) to the circuit

        This method takes care that its nodes are added as well.

        **Parameters:**

        part_id : string
            The cccs ID (eg ``'H1'``). The first letter is always ``'H'``.
        n1, n2 : strings
            The output port nodes, where the output current is
            forced. Eg. "outp", "outm" or "out_a", "out_b".
        source_id : string
            The voltage source to be used to sense the current that drives
            the output voltage. Eg. ``'V1'``.
        value : float
            The proportionality factor between the sense current :math:`I_s`
            flowing into the ``source_id`` voltage source (input) and output voltage.
            Mathematically:

            .. math::

                Vn_1 - Vn_2 = \\alpha I_s

        .. seealso::

            :class:`ahkab.components.sources.EVSource`,
            :class:`ahkab.components.sources.FISource`

        """
        # Add the nodes, this is SAFE: if a node is already known to the circuit,
        # the methods will just ignore the request.
        n1 = self.add_node(n1)
        n2 = self.add_node(n2)
        # instantiate the element
        elem = components.sources.HVSource(part_id=part_id, n1=n1, n2=n2,
                                source_id=source_id, value=value)
        # add it!
        self.append(elem)

    def add_vcvs(self, part_id, n1, n2, sn1, sn2, value):
        """Adds a voltage-controlled voltage source (vcvs) to the circuit

        This method also takes care that its nodes are added as well.

        **Parameters:**

        part_id : string
            The vcvs ID (eg "E1"). The first letter is always E.
        n1, n2 : string
            The output port nodes, where the output voltage is
            forced. Eg. "outp", "outm" or "out_a", "out_b".
        sn1, sn2 : string
            The input port nodes, where the input voltage is
            read. Eg. "inp", "inm" or "in_a", "in_b".
        alpha : float
            The proportionality factor between input and output voltages is
            given by the relationship:

            .. math::

                V(out_p) - V(out_n) = \\alpha \\cdot (V(in_p) - V(in_n))

        """

        n1 = self.add_node(n1)
        n2 = self.add_node(n2)
        sn1 = self.add_node(sn1)
        sn2 = self.add_node(sn2)

        elem = components.sources.EVSource(part_id=part_id, n1=n1, n2=n2, sn1=sn1, sn2=sn2,
                                value=value)

        self.append(elem)


    def add_vccs(self, part_id, n1, n2, sn1, sn2, value):
        """Adds a voltage-controlled current source (VCCS) to the circuit

        This method also takes care that its nodes are added as well.

        **Parameters:**

        part_id : string
            The VCCS ID (eg ``"G1"``). The first letter is always ``'G'``.
        n1, n2 : string
            The output port nodes, where the output current is
            forced. Eg. "outp", "outm" or "out_a", "out_b".
            The passive convention is used as everywhere else in the simulator:
            a positive current flows into ``n1`` and out of ``n2``.
        sn1, sn2 : string
            The input port nodes, where the input voltage is
            sensed. Eg. "inp", "inm" or "in_a", "in_b".
        value : float
            The proportionality factor between input and output voltages,
            which are related by the equality:

            .. math::

                I_o = alpha * \\left[V(inp) - V(inn)\\right]

        """

        n1 = self.add_node(n1)
        n2 = self.add_node(n2)
        sn1 = self.add_node(sn1)
        sn2 = self.add_node(sn2)

        elem = components.sources.GISource(part_id=part_id, n1=n1, n2=n2, sn1=sn1, sn2=sn2,
                                value=value)

        self.append(elem)

    def add_switch(self, part_id, n1, n2, sn1, sn2, ic, model_label, models=None):
        """Adds a voltage-controlled or current-controlled switch to the circuit

        This method also takes care that its nodes are added to the circuit as
        well, if necessary.

        **Notice:**

        - Current-controlled switches are not yet implemented. If you try to add
          one, you'll trigger an error. If you got a bit of time to spare,
          patches are welcome.
        - The switches ``part_id`` should begin with ``'S'`` for
          voltage-controlled switches and with ``'W'`` for current-controlled
          switches.
        - The actual behavior is set by the model. Make sure you supply a
          voltage-controlled switch model for a voltage-controlled switch and
          the appropriate type of model for the current-controlled switch.
          Mixing them up will go undetected.

        **Parameters:**

        part_id : string
            the switch ID (eg ``"S1"`` - voltage-controlled - or ``"Wa"`` -
            current-controlled). The first letter is always ``S`` or ``W``.
        n1, n2 : string
            the output port nodes, where the switch is connected. Eg. ``"out1"``,
            ``"out2"`` or ``"n_a"``, ``"n_b"``.
        sn1, sn2 : string
            The input port nodes, where the input voltage is
            read. Eg. "inp", "inm" or "in_a", "in_b".
        ic : boolean
            The initial conditions for transient simulation. Not currently
            implemented!
        model_label : string
            The switch model identifier. The model needs to be added
            first, then the elements using it.
        models : dict, optional
            A dictionary assembled as (identifier:instance), containing all the available model
            instances. If not set or ``None``, the circuit models will be used (recommended).
        """

        n1 = self.add_node(n1)
        n2 = self.add_node(n2)
        sn1 = self.add_node(sn1)
        sn2 = self.add_node(sn2)

        if models is None:
            models = self.models
        if model_label not in models:
            raise ModelError("Unknown switch model id: " + model_label)

        elem = switch.switch_device(part_id=part_id, n1=n1, n2=n2, sn1=sn1,
                                    sn2=sn2, model=models[model_label])

        self.append(elem)

    def add_user_defined(self, module, label, param_dict):
        """Adds a user defined element.

        In order for this to work, you should write a module that supplies the
        elem class.

        XXX WRITE DOC
        """

        if module_name in circuit.user_defined_modules_dict:
            module = circuit.user_defined_modules_dict[module_name]
        else:
            fp, pathname, description = importlib.util.find_spec(module_name)
            module = importlib.import_module(module_name, fp, pathname, description)
            circuit.user_defined_modules_dict.update({module_name: module})

        elem_class = getattr(module, label)

        param_dict.update({"convert_units": convert_units})
        param_dict.update({"circuit_node": self.add_node})

        elem = elem_class(**param_dict)
        elem.part_id = "y%s" % part_id[1:]

        # call check() if supported
        if hasattr(elem, "check"):
            selfcheck_result, error_msg = elem.check()
            if not selfcheck_result:
                raise NetlistParseError("module: " + module_name + \
                                        " elem type: " + elem_type_name + \
                                        " error: " + error_msg)

        self.append(elem)

    def remove_elem(self, elem_or_id):
        """Removes an element from the circuit and takes care that no
        "orphan" nodes are left.

        .. note::

            Support for removing elements is experimental.

        **Parameters:**

        elem_or_id : string or circuit element
            You may pass as first element, alternatively, either the ``part_id``
            of the element to be removed or the element itself.

        The method will also take care of purging from the circuit nodes that
        are left orphan, ie with no elements connected.

        :raises ValueError: if no such element is found in the circuit.

        """
        if type(elem_or_id) in py3compat.string_types:
            #we got a part_id, we need the element
            elem = self.get_elem_by_name(elem_or_id)
        else:
            # we got the element
            elem = elem_or_id
        self.remove(elem)
        nodes = []
        if hasattr(elem, 'n1') and elem.n1 != 0:
            nodes = nodes + [elem.n1]
        if hasattr(elem, 'n2') and elem.n2 != 0 and elem.n2 not in nodes:
            nodes = nodes + [elem.n2]
        if elem.is_nonlinear:
            for n1, n2 in elem.ports:
                if n1 != 0 and n1 not in nodes:
                    nodes = nodes + [n1]
                if n2 != 0 and n2 not in nodes:
                    nodes = nodes + [n2]
        # check if then nodes are needed by other elements
        remove_list = copy.copy(nodes)
        for n in nodes:
            for e in self:
                if hasattr(elem, 'n1') and e.n1 == n or \
                   hasattr(elem, 'n2') and e.n2 == n or \
                   hasattr(elem, 'sn1') and e.sn1 == n or \
                   hasattr(elem, 'sn2') and e.sn2 == n:
                    remove_list.remove(n)
                    break
                if elem.is_nonlinear:
                    oports = elem.get_output_ports()
                    # check output ports
                    for n1, n2 in oports:
                        if n1 == n or n2 == n:
                            remove_list.remove(n)
                            break
                    if n not in remove_list:
                        break
                    # check the ports that drive them
                    for i in range(len(oports)):
                        dports = elem.get_drive_ports(i)
                        for n1, n2 in oports:
                            if n1 == n or n2 == n:
                                remove_list.remove(n)
                                break
                        if n not in remove_list:
                            break

        for n in remove_list:
            self.nodes_dict.pop(self.nodes_dict[n])
            self.nodes_dict.pop(n)

    def find_vde_index(self, elem_or_id, verbose=3):
        """Finds a voltage-defined element MNA index.

        **Parameters:**

        elem_or_id : string or circuit element
            You may pass as first element, alternatively, either the ``part_id``
            of the element whose index is being requested (eg. 'V1') or the
            element itself.
            Notice the ``part_id`` includes both the id letter (eg. 'V') and the
            description (eg. '1').
        verbose : int
            The verbosity level, from 0 (silent) to 6 (debug).

        **Returns:**

        indx : int
            The index.

        :raises ValueError: if no such element is in the circuit.
        """
        if type(elem_or_id) not in py3compat.string_types:
            # we got an element
            part_id = elem_or_id.part_id
        else:
            # we got a string corresponding to the part_id of an element
            part_id = elem_or_id
        vde_index = 0
        for elem in self:
            if is_elem_voltage_defined(elem):
                if elem.part_id.upper() == part_id.upper():
                    break
                else:
                    vde_index += 1
        else:
            raise ValueError(("find_vde_index(): element %s was not found." +\
                              " This is a bug.") % (part_id,))
        printing.print_info_line(("%s found at index %d" % (part_id,
                                                            vde_index), 6),
                                 verbose)
        return vde_index

    def find_vde(self, index):
        """Finds a voltage-defined element from its MNA KVL index

        **Parameters:**

        index : int
            The element index in the KVL equations.

        **Returns:**

        e : circuit element (an instance of a subclass of Component)
            The element corresponding to ``index``.

        :raises IndexError: if no element corresponds to such an index.
        """
        index = index - len(self.nodes_dict)/2 + 1
        ni = 0
        for e in self:
            if is_elem_voltage_defined(e):
                if index == ni:
                    break
                else:
                    ni = ni + 1
        else: #executed if no break occurred
            raise IndexError('No element corresponds to vde index %d' %
                             (index + len(self.nodes_dict)/2 - 1))
        return e


# STATIC METHODS
def is_elem_voltage_defined(elem):
    """Check if an element needs its own KCL equation

    **Parameters:**

    elem : Component
        The element to be checked.

    **Returns:**

    chk : bool
        ``True`` if ``elem`` is a voltage source, an inductor, a voltage-controlled
        voltage source or a current-controlled voltage source. ``False`` otherwise.
    """
    if isinstance(elem, components.sources.VSource) or isinstance(elem, components.sources.EVSource) or \
        isinstance(elem, components.sources.HVSource) or isinstance(elem, components.Inductor) \
            or (hasattr(elem, "is_voltage_defined") and elem.is_voltage_defined):
        return True
    else:
        return False


class NodeNotFoundError(Exception):
    """Circuit Node exception."""
    pass


class CircuitError(Exception):
    """General circuit assembly exception."""
    pass


class ModelError(Exception):
    """Model not found exception."""
    pass


class subckt:
    """This class holds the necessary information about a circuit.

    An instance of this class is returned by:

    :func:`ahkab.netlist_parser.parse_sub_declaration`

    **Parameters:**

    name : string
        The subcircuit definition label.
    code : string
        The netlist code that can be instantiated have a circuit
        instance.
    connected_nodes_list : list
        A list of nodes that are used in the circuit and that are
        meant to be connected to the external circuit.

    Notice that in the current implementation, the GND node (0)
    is *always* global.

    """
    def __init__(self, name, code, connected_nodes_list):
        self.name = name
        self.connected_nodes_list = connected_nodes_list
        self.code = code


class _circuit_wrapper:
    """Fictious circuit class, meant to wrap subcircuits.

    Not meant for end users at this stage.

    Rationale:

    Within a subcircuit, the nodes name are fictious.
    All nodes have to be renamed before a subcircuit instance
    may be insterted in a circuit (all our circuits are flat in
    memory for now).

    The nodes of the subcircuit that are connected to the
    nodes of the circuit have to be renamed to them, those
    that are not referenced there need to be renamed in order
    for them ot be unique.

    This class wraps a circuit object and performs the conversion
    _before_ calling ``circ.add_node()``.

    While instatiating/calling a subcircuit wrap the circuit instance
    in this.

    **Parameters:**

    circ : circuit instance
        The main circuit. Remember that all our assembled circuits
        are flat in memory.
    connection_nodes_dict : dictionary
        The dictionary mapping internal nodes to global, circuit-wide nodes.
    subckt_name : string
        The subcircuit instance name. The first letter must always be ``'X'``.
    subckt_label : string
        The label of the subcircuit that is being instantiated.

    """

    def __init__(self, circ, connection_nodes_dict, subckt_name, subckt_label):
        self.circ = circ
        self.prefix = subckt_label + "-" + subckt_name + "-"
        self.subckt_node_filter_dict = {}
        self.subckt_node_filter_dict.update(connection_nodes_dict)
        self.subckt_node_filter_dict.update({'0': '0'})

    def add_node(self, ext_name):
        """We want to perform the following conversions:

        * connected node in the subcircuit -> node in the upper circuit
        * local-only node of the subcircuit -> rename it to something unique
        * REF (0) -> REF (0)

        And then call ``circ.add_node()``.
        """
        if ext_name not in self.subckt_node_filter_dict:
            self.subckt_node_filter_dict.update(
                {ext_name: self.prefix + ext_name})
        int_node = self.circ.add_node(self.subckt_node_filter_dict[ext_name])
        return int_node
