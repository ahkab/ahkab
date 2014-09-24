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

Every circuit is described in the ahkab simulator by an instance of the
Circuit class.
This class holds everything is needed to simulate the circuit (except
the specification of the analyses to be performed).

It is even possible to rewrite a netlist from a Circuit class: see the
printing module.

The Circuit
\"\"\"\"\"\"\"\"\"\"\"

A circuit is derived from a list. This list holds the elements.

All the elements in the ``Circuit`` class must be appended to their
``Circuit`` instance and connections should be ensure checking that
the nodes the elements refer to are indeed circuit nodes..

To simplify the operation of adding a component to a ``Circuit``,
the following convenience methods are provided to the user to add and
remove elements to the circuit:

* :func:`Circuit.add_resistor`

* :func:`Circuit.add_capacitor`

* :func:`Circuit.add_inductor`

* :func:`Circuit.add_vsource`

* :func:`Circuit.add_isource`

* :func:`Circuit.add_diode`

* :func:`Circuit.add_mos`

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
    mycircuit.add_resistor(name="R1", n1="n1", n2=gnd, R=600)

Refer to the methods help for addtional info.

Nodes
\"\"\"\"\"

The nodes are stored in this way: we assign to each node a internal
name, whatever is its external one (which is used in the netlist).
Those are integers.

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

They are stored in ``Circuit.models`` (type dict), the following methods
are provided to add and remove device models.

* :func:`Circuit.add_model`

* :func:`Circuit.remove_model`

Reference
\"\"\"\"\"\"\"\"\"

"""

import sys
import math

from . import devices
from . import diode
from . import ekv
from . import mosq
from . import printing
from . import switch

# will be added here by netlist_parser and circuit instances
user_defined_modules_dict = {}

class Circuit(list):
    """The circuit class.

    **Parameters:**

    title : string
        The circuit title.

    filename : string, optional
        If the circuit instance corresponds to a netlist file on disk,
        set this to the netlist filename.

    .. deprecated:: the filename option was deprecated in v. 0.09
    """
    def __init__(self, title, filename=None):
        self.title = title
        self.filename = filename
        self.nodes_dict = {}  # {int_node:ext_node}
        self.internal_nodes = 0
        self.models = {}
        self.gnd = '0'

    def __str__(self):
        s = "* " + self.title + "\n"
        for elem in self:
            s += elem.get_netlist_elem_line(self.nodes_dict) + "\n"
        return s[:-1]

    def create_node(self, name):
        """Creates a new node, adds it to the circuit and returns it to the user
        (to be used for subsequent declaration of elements, for example).

        If there is a node in the circuit with the same name, ValueError is
        raised.

        **Parameters:**

        name : string
            the _unique_ identifier of the node.

        **Returns:**

        node : string 
            the _unique_ identifier of the node.

        :raises ValueError: if a new node with the given id cannot be created,
          for example because a node with the same name already exists in the
          circuit. The only exception is the ground node, which has the 
          reserved id ``'0'``, and for which this method won't raise any
          exception.

        """
        got_ref = self.nodes_dict.has_key(0)
        try:
            self.nodes_dict.values().index(name)
        except ValueError:
            if name == '0':
                int_node = 0
            else:
                int_node = len(self.nodes_dict) + 1 * (not got_ref)
            self.nodes_dict.update({int_node: name})
        else:
            raise ValueError
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
        for the circuit to be non-patological and has ``ext_name=str(int_name)='0'``.

        Notice that this method doesn't halt or print errors if the node is already been
        added previsiously. It simply returns the internal node name assigned to it.

        **Parameters:**

        ext_name : string
            The unique identifier of the node.

        **Returns:**

        int_node : int
            the internal node id assigned to the node.

        """
        got_ref = self.nodes_dict.has_key(0)

        # test: do we already have it in the dictionary?
        try:
            self.nodes_dict.values().index(ext_name)
        except ValueError:
            if ext_name == '0':
                int_node = 0
            else:
                int_node = len(self.nodes_dict) + 1 * (not got_ref)
            self.nodes_dict.update({int_node: ext_name})
        else:
            for (key, value) in self.nodes_dict.iteritems():
                if value == ext_name:
                    int_node = key
        return int_node

    def generate_internal_only_node_label(self):
        """Generate implicit internal nodes.

        Some devices are made of a group of other devices, connected by
        "internal only" nodes, which have the prefix ``'INT'`` and the
        simulator treats specially, hiding them from the user if not
        explicitly asked about them.

        This method generates the external names for such nodes. 

        .. note:: 

            They are *NOT* added to the circuit during their creation.

        **Returns:**

        ext_node : string
            The corresponding external node id.
        """

        ext_node = "INT" + str(self.internal_nodes)
        self.internal_nodes = self.internal_nodes + 1
        return ext_node

    def is_int_node_internal_only(self, int_node):
        """Check whether a supplied node id corresponds to an "internal node"
        or not.
        """
        return self.nodes_dict[int_node].find("INT") > -1

    def is_nonlinear(self):
        """Check whether the circuit is non-linear or not.
        """
        for elem in self:
            if elem.is_nonlinear:
                return True
        return False

    def get_locked_nodes(self):
        """Get all nodes connected to non-linear elements.

        This list is meant to be passed to dc_solve or mdn_solver to be used in
        get_td to evaluate the damping coefficient in a Newton-Rhapson
        iteration.

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

        .. note::

            This method is slow, because it has to look through ``Circuit.nodes_dict``.

        Throws a NodeNotFoundError exception.


        **Returns:**

        int_node : int
            The internal node associated.

        """
        items = self.nodes_dict.items()
        values = [value for key, value in items]

        try:
            index = values.index(ext_node)
        except ValueError:
            raise NodeNotFoundError, "Node %s not found in the circuit." % ext_node

        return items[index][0]

    def int_node_to_ext(self, int_node):
        """This function returns the string id associated with the integer internal node id
        ``int_node``.

        **Parameters:**

        int_node : int
            The internal node id to be converted.

        .. note:: 

            Accessing this function is the same as referencing ``circuit_inst.nodes_dict[int_node]``, 
            except a ``NodeNotFoundError`` exception is thrown instead of a ``KeyError``.

        .. note::

            This method is fast much faster than :func:`Circuit.ext_node_to_int`.

        **Returns:** 

        ext_node : string 
            the string id associated with ``int_node``.
        """
        try:
            ret = self.nodes_dict[int_node]
        except KeyError:
            raise NodeNotFoundError, ""

        return ret

    def has_duplicate_elem(self):
        """Self check.

        No circuit should ever have duplicate elements
        (ie elements with the same ``part_id``)."""
        for index1 in range(len(self)):
            for index2 in range(index1 + 1, len(self)):
                if self[index1].part_id == self[index2].part_id:
                    return True
        return False

    def get_ground_node(self):
        """Returns the reference node, AKA gnd."""
        return '0'

    def get_elem_by_name(self, name):
        """Get a circuit element from its ``part_id`` value.

        If no matching element is found, the method returns
        ``None``. This may change in the future.

        **Parameters:**

        name : string
            The ``part_id`` of the element 

        **Returns:**

        elem : circuit element or None
            Depending whether a matching element was found or not.
        """
        for e in self:
            if e.part_id.lower() == name.lower():
                return e
        return None

    def add_model(self, model_type, model_label, model_parameters):
        """Add a model to the available circuit models.

        **Parameters:**

        model_type : string
            the model type (eg "ekv").

        model_label : string
            a unique identifier for the model being added (eg. "nch1").

        model_parameters: dict
            a dictionary holding the parameters to be supplied to the 
            model to initialize it.

        **Returns:**

        models : list
            the updated models

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
            raise CircuitError, "Unknown model type %s" % (model_type,)
        self.models.update({model_label: model_iter})
        return self.models

    def remove_model(self, model_label):
        """Remove a model from the available models.

        **Parameters:**

        model_label : string
            the unique identifier corresponding to the model
            being removed.

        .. note::

            This method currently silently ignores models that are not defined.

        """
        if self.models is not None and self.models.has_key(model_label):
            del self.models[model_label]
        # should print a warning here

    def add_resistor(self, name, n1, n2, value):
        """Adds a resistor to the circuit.

        The resistor instance is added to the circuit elements 
        and connected to the provided nodes. If the nodes are not
        found in the circuit, they are created and added as well.

        **Parameters:**

        name : string
            the resistor name (eg "R1"). The first letter is replaced by an R

        n1, n2 : string
            the nodes to which the resistor is connected.

        value : float,
            The resistance between ``n1`` and ``n2`` in Ohm.

        **Returns:** 

            True

        .. seealso:: 

            :func:`add_resistor`, :func:`add_capacitor`, 
            :func:`add_inductor`, :func:`add_vsource`, :func:`add_isource`,
            :func:`add_diode`, :func:`add_mos`, :func:`add_vcvs`, :func:`add_vccs`,
            :func:`add_user_defined`, :func:`remove_elem`
 
        """
        n1 = self.add_node(n1)
        n2 = self.add_node(n2)

        if value == 0:
            raise CircuitError, "ZERO-valued resistors are not allowed."

        elem = devices.Resistor(n1=n1, n2=n2, value=value)
        elem.part_id = name
        self.append(elem)
        return True

    def add_capacitor(self, name, n1, n2, value, ic=None):
        """Adds a capacitor to the circuit. 

        The capacitor instance is added to the circuit elements 
        and connected to the provided nodes. If the nodes are not
        found in the circuit, they are created and added as well.

        **Parameters:**

        name : string
            The capacitor name (eg "C1"). The first letter is always C.

        n1, n2 : string
            The nodes to which the element is connected.

        value : float
            The capacitance value.

        ic : float, optional
            The initial condition, if any. See the simulation docs for
            how this affects the results.

        Returns: True

        .. seealso:: 

        :func:`add_resistor`,
        :func:`add_inductor`, :func:`add_vsource`, :func:`add_isource`,
        :func:`add_diode`, :func:`add_mos`, :func:`add_vcvs`, :func:`add_vccs`,
        :func:`add_user_defined`, :func:`remove_elem`

        """
        if value == 0:
            raise CircuitError, "ZERO-valued capacitors are not allowed."

        n1 = self.add_node(n1)
        n2 = self.add_node(n2)

        elem = devices.Capacitor(n1=n1, n2=n2, value=value, ic=ic)
        elem.part_id = name

        self.append(elem)
        return True

    def add_inductor(self, name, n1, n2, value, ic=None):
        """Adds an inductor to the circuit. 

        The inductor instance is added to the circuit elements 
        and connected to the provided nodes. If the nodes are not
        found in the circuit, they are created and added as well.

        **Parameters:**

        name : string
            the inductor name (eg "Lfilter"). The first letter is always L.

        n1, n2 : string
            the nodes to which the element is connected. Eg. ``"in"`` or ``"out_a"``.

        value : float
            The inductance value.

        ic : float, optional
            Initial condition, see simulation types for how this affects
            the results.

        **Returns:** 

        True

        .. seealso:: 

            :func:`add_resistor`, :func:`add_capacitor`, 
            :func:`add_inductor`, :func:`add_vsource`, :func:`add_isource`,
            :func:`add_diode`, :func:`add_mos`, :func:`add_vcvs`, :func:`add_vccs`,
            :func:`add_user_defined`, :func:`remove_elem`
        """

        n1 = self.add_node(n1)
        n2 = self.add_node(n2)

        elem = devices.Inductor(n1=n1, n2=n2, value=value, ic=ic)
        elem.part_id = name

        self.append(elem)
        return True

    def add_inductor_coupling(self, name, L1, L2, value):
        """ Write DOC XXX
        """
        L1elem, L2elem = None, None

        for e in self:
            if isinstance(e, devices.Inductor) and (L1 == e.part_id):
                L1elem = e
            elif isinstance(e, devices.Inductor) and (L2 == e.part_id):
                L2elem = e

        if L1elem is None or L2elem is None:
            error_msg = "One or more coupled inductors for %s were not found: %s (found: %s), %s (found: %s)." % \
                (name, L1, L1elem is not None, L2, L2elem is not None)
            printing.print_general_error(error_msg)
            printing.print_general_error("Quitting.")
            sys.exit(30)

        M = math.sqrt(L1elem.value * L2elem.value) * value

        elem = devices.InductorCoupling(L1=L1, L2=L2, K=value, M=M)
        elem.part_id = name
        L1elem.coupling_devices.append(elem)
        L2elem.coupling_devices.append(elem)

        self.append(elem)

    def add_vsource(self, part_id, n1, n2, dc_value, ac_value=0, function=None):
        """Adds a voltage source to the circuit (also takes care that the nodes
        are added as well).

        Parameters:
        name (string): the voltage source name (eg "VA"). The first letter is always V.
        n1, n2 (string): the nodes to which the element is connected.
                    eg. "in" or "out_a"
        dc_value (float): DC voltage
        ac_value (float): AC voltage (optional)
        function (function): optional time function. See devices.py for built-ins.

        Returns: True
        """
        n1 = self.add_node(n1)
        n2 = self.add_node(n2)

        elem = devices.VSource(
            part_id=part_id, n1=n1, n2=n2, dc_value=dc_value, ac_value=ac_value)

        if function is not None:
            elem.is_timedependent = True
            elem._time_function = function

        self.append(elem)
        return True

    def add_isource(self, part_id, n1, n2, dc_value, ac_value=0, function=None):
        """Adds a current source to the circuit (also takes care that the nodes
        are added as well).

        Parameters:
        name (string): the current source name (eg "IA"). The first letter is always I.
        n1, n2 (string): the nodes to which the element is connected.
                    eg. "in" or "out_a"
        dc_value (float): DC current
        ac_value (float): AC current (optional)
        function (function): optional time function. See devices.py for built-ins.

        Returns: True
        """
        n1 = self.add_node(n1)
        n2 = self.add_node(n2)

        elem = devices.ISource(
            part_id=part_id, n1=n1, n2=n2, dc_value=dc_value, ac_value=ac_value)

        if function is not None:
            elem.is_timedependent = True
            elem._time_function = function

        self.append(elem)
        return True

    def add_diode(self, part_id, n1, n2, model_label, models=None, Area=None, T=None, ic=None, off=False):
        """Adds a diode to the circuit (also takes care that the nodes
        are added as well).

        Parameters:
        name (string): the diode name (eg "D1"). The first letter is always D.
        n1, n2 (string): the nodes to which the element is connected.
                    eg. "in" or "out_a"
        Area (float): Scaled device area (optional, defaults to 1.0)
        T (float): operating temperature (no temperature dependence yet)
        ic (float): initial condition (not really implemented yet)
        model_label (string): the diode model identifier. The model needs to be added
                              first, then the elements using it.
        models (dict(identifier:instance), optional): list of available model
            instances. If not set or None, the circuit models will be used (recommended).

        Returns: True
        """
        n1 = self.add_node(n1)
        n2 = self.add_node(n2)
        if models is None:
            models = self.models
        if not models.has_key(model_label):
            raise ModelError, "Unknown diode model id: " + model_label

        elem = diode.diode(part_id=part_id, n1=n1, n2=n2, model=models[
                           model_label], AREA=Area, T=T, ic=ic, off=off)
        self.append(elem)

        return True

    def add_mos(self, part_id, nd, ng, ns, nb, w, l, model_label, models=None, m=None, n=None):
        """Adds a mosfet to the circuit (also takes care that the nodes
        are added as well).

        Parameters:

        name (string): the mos name (eg "M1"). The first letter is always M.
        nd (string): drain node
        ng (string): gate node
        ns (string): source node
        nb (string): bulk node
        w (float): gate width
        l (float): gate length
        model_label (string): model identifier
        models (circuit models): circuit models
        m (int): shunt multiplier (optional)
        n (int): series multiplier (unsupported)

        Returns: True
        """
        if m is None:
            m = 1
        if n is None:
            n = 1

        nd = self.add_node(nd)
        ng = self.add_node(ng)
        ns = self.add_node(ns)
        nb = self.add_node(nb)

        if models is None:
            models = self.models

        if not models.has_key(model_label):
            raise ModelError, "Unknown model id: " + model_label

        if isinstance(models[model_label], ekv.ekv_mos_model):
            elem = ekv.ekv_device(
                part_id, nd, ng, ns, nb, w, l, models[model_label], m, n)

        elif isinstance(models[model_label], mosq.mosq_mos_model):
            elem = mosq.mosq_device(
                nd, ng, ns, nb, w, l, models[model_label], m, n, part_id)

        else:
            raise Exception, "Unknown model type for " + model_label

        self.append(elem)

        return True

    def add_vcvs(self, part_id, n1, n2, sn1, sn2, value):
        """Adds a voltage-controlled voltage source (vcvs) to the circuit
        (also takes care that its nodes are added as well).

        Parameters:
        name (string): the vcvs name (eg "E1"). The first letter is always E.
        n1, n2 (string): the output port nodes, where the output voltage is
                     forced. Eg. "outp", "outm" or "out_a", "out_b".
        sn1, sn2 (string): the input port nodes, where the input voltage is
                     read. Eg. "inp", "inm" or "in_a", "in_b".
                alpha (float): The proportionality factor between input and output voltages:
                V(outp) - V(outn) = alpha * (V(inp) - V(inn))

        Returns: True
        """

        n1 = self.add_node(n1)
        n2 = self.add_node(n2)
        sn1 = self.add_node(sn1)
        sn2 = self.add_node(sn2)

        elem = devices.EVSource(
            part_id=part_id, n1=n1, n2=n2, sn1=sn1, sn2=sn2, value=value)

        self.append(elem)

        return True

    def add_vccs(self, part_id, n1, n2, sn1, sn2, value):
        """Adds a voltage-controlled current source (vccs) to the circuit
        (also takes care that its nodes are added as well).

        Parameters:
        name (string): the vccs name (eg "G1"). The first letter is always G.
        n1, n2 (string): the output port nodes, where the output current is
                     forced. Eg. "outp", "outm" or "out_a", "out_b".
                     The usual convention is used: a positive current
                     flows into n1 and out of n2.
        sn1, sn2 (string): the input port nodes, where the input voltage is
                       read. Eg. "inp", "inm" or "in_a", "in_b".
                alpha (float): The proportionality factor between input and output voltages:
                I[G1] = alpha * (V(inp) - V(inn))

        Returns: True
        """

        n1 = self.add_node(n1)
        n2 = self.add_node(n2)
        sn1 = self.add_node(sn1)
        sn2 = self.add_node(sn2)

        elem = devices.GISource(
            part_id=part_id, n1=n1, n2=n2, sn1=sn1, sn2=sn2, value=value)

        self.append(elem)
        return True

    def add_switch(self, name, n1, n2, sn1, sn2, ic, model_label, models=None):
        """Adds a voltage-controlled or current-controlled switch to the circuit
        (also takes care that its nodes are added as well).

        Notice:
        - Current-controlled switches are not yet implemented. If you try to add one
          you'll trigger an error.
        - The switches name should begin with 'S' for voltage-controlled switches
          and with 'W' for current-controlled switches.
        - The actual behavior is set by the model. Make sure you supply a voltage-controlled
          switch model for a voltage-controlled switch and the same for the
          current-controlled switch. Mixing them up will go undetected.

        Parameters:
        name (string): the switch name (eg "S1" - voltage-controlled - or 'W1' -
                       current-controlled). The first letter is always S or W.
        n1, n2 (string): the output port nodes, where the switch is connected.
                     Eg. "outp", "outm" or "out_a", "out_b".
        sn1, sn2 (string): the input port nodes, where the input voltage is
                       read. Eg. "inp", "inm" or "in_a", "in_b".
        ic (boolean): the initial conditions for transient simulation. Not currently
                      implemented!
        model_label (string): the switch model identifier. The model needs to be added
                              first, then the elements using it.
        models (dict(identifier:instance), optional): list of available model
            instances. If not set or None, the circuit models will be used (recommended).

        Returns: True
        """

        n1 = self.add_node(n1)
        n2 = self.add_node(n2)
        sn1 = self.add_node(sn1)
        sn2 = self.add_node(sn2)

        if models is None:
            models = self.models
        if not models.has_key(model_label):
            raise ModelError, "Unknown switch model id: " + model_label

        elem = switch.switch_device(
            part_id=part_id, n1=n1, n2=n2, sn1=sn1, sn2=sn2, model=models[model_label])
        self.append(elem)
        return True

    def add_user_defined(self, module, label, param_dict):
        """Adds a user defined element.

        In order for this to work, you should write a module that supplies the
        elem class.

        XXX WRITE DOC
        """

        if circuit.user_defined_modules_dict.has_key(module_name):
            module = circuit.user_defined_modules_dict[module_name]
        else:
            fp, pathname, description = imp.find_module(module_name)
            module = imp.load_module(module_name, fp, pathname, description)
            circuit.user_defined_modules_dict.update({module_name: module})

        elem_class = getattr(module, label)

        param_dict.update({"convert_units": convert_units})
        param_dict.update({"circuit_node": self.add_node})

        elem = elem_class(**param_dict)
        elem.part_id = "y%s" % name[1:]

        if hasattr(elem, "check"):
            selfcheck_result, error_msg = elem.check()
            if not selfcheck_result:
                raise NetlistParseError, "module: " + module_name + " elem type: " + elem_type_name + " error: " +\
                    error_msg

        self.append(elem)
        return True

    def remove_elem(self, elem):
        """Removes an element from the circuit and takes care that no
        "orphan" nodes are left.
        circ: the circuit instance
        elem: the element to be removed

        Returns: True if the element was found and removed, False otherwise
        """
        if elem not in self:
            return False

        self.remove(elem)
        nodes = []
        if hasattr(elem, n1) and elem.n1 != 0:
            nodes = nodes + [n1]
        if hasattr(elem, n2) and elem.n2 != 0 and elem.n2 not in nodes:
            nodes = nodes + [n2]
        if elem.is_nonlinear:
            for n1, n2 in elem.ports:
                if n1 != 0 and n1 not in nodes:
                    nodes = nodes + [n1]
                if n2 != 0 and n2 not in nodes:
                    nodes = nodes + [n2]

        remove_list = copy.copy(nodes)
        for n in nodes:
            for e in self:
                if hasattr(elem, n1) and e.n1 == n or\
                        hasattr(elem, n2) and e.n2 == n:
                    remove_list.remove(n)
                    break
                if elem.is_nonlinear:
                    for n1, n2 in elem.ports:
                        if n1 == n or n2 == n:
                            remove_list.remove(n)
        for n in remove_list:
            self.nodes_dict.pop(n)
        return True

    def find_vde_index(self, id_wdescr, verbose=3):
        """Finds a voltage defined element MNA index.

        Parameters:
        id_wdescr (string): the element name, eg. 'V1'. Notice it includes
                            both the id ('V') and the description ('1').
        verbose (int): verbosity level, from 0 (silent) to 6 (debug).

        Returns:
        the index (int)
        """
        vde_index = 0
        found = False
        for elem in self:
            if is_elem_voltage_defined(elem):
                if elem.part_id.upper() == id_wdescr.upper():
                    found = True
                    break
                else:
                    vde_index += 1

        if not found:
            printing.print_warning(
                "find_vde_index(): element %s was not found. This is a bug." % (id_wdescr,))
        else:
            printing.print_info_line(
                ("%s found at index %d" % (id_wdescr, vde_index), 6), verbose)
        return vde_index

    def find_vde(self, index):
        """Finds a voltage defined element MNA index.

        Parameters:
        id_wdescr (string): the element name, eg. 'V1'. Notice it includes
                            both the id ('V') and the description ('1').
                            The search term is case insensitive.

        Returns:
        the index (int)
        """
        index = index - len(self.nodes) + 1
        ni = 0
        found = False
        for e in self:
            if circuit.is_elem_voltage_defined(e):
                if index == ni:
                    found = True
                else:
                    ni = ni + 1
                if found:
                    break
        if found:
            ret = e
        else:
            ret = None
        return ret


# STATIC METHODS
def is_elem_voltage_defined(elem):
    """Returns:
    True if the elem is a vsource, inductor, evsource or hvsource
    False otherwise.
    """
    if isinstance(elem, devices.VSource) or isinstance(elem, devices.EVSource) or \
        isinstance(elem, devices.HVSource) or isinstance(elem, devices.Inductor) \
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

      netlist_parser.parse_sub_declaration(subckt_lines)


    """
    name = ""
    connected_nodes_list = []
    code = []

    def __init__(self, name, code, connected_nodes_list):
        self.name = name
        self.connected_nodes_list = connected_nodes_list
        self.code = code


class circuit_wrapper:

    """Within a subcircuit, the nodes name are fictious.
    The nodes of the subcircuit that are connected to the
    nodes of the circuit have to be renamed to them, the
    others have to be renamed too.

    This class wraps a circuit object and performs the conversion
    _before_ calling circ.add_node()

    While instatiating/calling a subcircuit wrap circ in this.
    """

    def __init__(self, circ, connection_nodes_dict, subckt_name, subckt_label):
        self.circ = circ
        self.prefix = subckt_label + "-" + subckt_name + "-"
        self.subckt_node_filter_dict = {}
        self.subckt_node_filter_dict.update(connection_nodes_dict)
        self.subckt_node_filter_dict.update({'0': '0'})

    def add_node(self, ext_name):
        """We want to perform the following conversions:
        connected node in the subcircuit -> node in the upper circuit
        local-only node of the subcircuit -> rename it to something uniq
        REF (0) -> REF (0)

        And then call circ.add_node()
        """
        if not self.subckt_node_filter_dict.has_key(ext_name):
            self.subckt_node_filter_dict.update(
                {ext_name: self.prefix + ext_name})
        int_node = self.circ.add_node(self.subckt_node_filter_dict[ext_name])
        return int_node
