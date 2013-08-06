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

import sys
import devices, printing, ekv, mosq

# will be added here by netlist_parser and circuit instances
user_defined_modules_dict = {}

class circuit:
	"""Every circuit is described in the ahkab simulator by a circuit class.
	This class holds everything is needed to simulate the circuit (except
	the specification of the analyses to be performed).
	
	It is even possible to rewrite a netlist from a circuit class: see the 
	printing module.

	There are basically three things in this class.
	
	1. Nodes
	The nodes are stored in this way: we assign to each node a internal 
	name, whatever is its external one (which is used in the netlist).
	Those are integers.

	The simulator uses always the internal names. When the results are
	presented to the user, the internal node is not showed, the external 
	identifier (or external node name) is printed instead.

	This is done through:
		my_circuit = circuit()
 
 		...
		[ init code ]
 		...

		print "This is a node" + my_circuit.nodes_dict[int_node]
	
	2. Elements
	All the elements in the circuit must be appended to the element list.

	The following methods are provided to add and remove elements to the circuit:

	add_resistor(self, name, ext_n1, ext_n2, R)
	add_capacitor(self, name, ext_n1, ext_n2, C, ic=None)
	add_inductor(self, name, ext_n1, ext_n2, L, ic=None)
	add_vsource(self, name, ext_n1, ext_n2, vdc, vac, function=None)
	add_isource(self, name, ext_n1, ext_n2, idc, iac, function=None)
	add_diode(self, name, ext_n1, ext_n2, Is=None, Rs=None, m=None, T=None, ic=None)
	add_mos(self, name, ext_nd, ext_ng, ext_ns, ext_nb, w, l, model_label, models, m=None, n=None)
	add_vcvs(self, name, ext_n1, ext_n2, ext_sn1, ext_sn2, alpha)
	add_vccs(self, name, ext_n1, ext_n2, ext_sn1, ext_sn2, alpha)
	add_user_defined(self, module, label, param_dict)
	remove_elem(self, elem)

	Example:
	
	mycircuit = circuit.circuit(title="Example circuit", filename=None)
	# no filename since there will be no deck associated with thios circuit.
	# get the ref node (gnd)
	gnd = mycircuit.get_ground_node()
	# add a node named n1 and a 600 ohm resistor connected between n1 and gnd
	mycircuit.add_resistor(name="R1", ext_n1="n1", ext_n2=gnd, R=600)
	
	Refer to the methods help for addtional info.

	3. Internal only nodes
	The number of internal only nodes (added automatically by the simulator)
	is held in my_circuit.internal_nodes
	That value shouldn't be changed by hand.

	4. Device models.
	They are stored in circuit.models (type dict), the following methods
	are provided to add and remove device models.

	add_model(self, model_type, model_label, model_parameters)
	remove_model(self, model_label)

	"""

	def __init__(self, title, filename):
		self.title = title
		self.filename = filename
		self.nodes_dict = {} # {int_node:ext_node}
		#_reverse_dict = {}
		self.elements = []
		self.internal_nodes = 0
		self.models = {}
	
	def add_node_to_circ(self, ext_name):
		"""Adds the supplied node to the circuit, if needed.
		
		When a 'normal' (not ref) node is added through add_node_to_circ(), a internal-only 
		name (or label) is assigned to it.
		
		The nodes labels are stored this way: self.nodes_dict is a dictionary of pairs 
		like (int_node:ext_node).
		
		Those internal names are integers, by definition, and are added starting from 1, 
		then 2,3,4,5...
		0 is reserved for the reference node (gnd), which is required and it has ext_name=='0'.
		
		Notice that this method doesn't halt or print errors if the node is already been
		added previsiously. It simply returns the internal node name assigned to it.
		
		Parameters:
		ext_name: a string that is used as _unique_ identifier of the node. 
		
		Returns: the internal node name (an INTEGER) assigned to the node.
		"""
		got_ref = self.nodes_dict.has_key(0)
		
		#test: do we already have it in the dictionary?
		try:
			self.nodes_dict.values().index(ext_name)
		except ValueError:
			if ext_name == '0':
				int_node = 0
			else:
				int_node = len(self.nodes_dict) + 1*(not got_ref)
			self.nodes_dict.update({int_node:ext_name})
		else:
			for (key, value) in self.nodes_dict.iteritems():
				if value == ext_name:
					int_node = key
		return int_node
		
	def generate_internal_only_node_label(self):
		"""Some devices are made of a group of other devices, connected by "internal only" nodes. 
		This method generates the external names for such nodes. They are NOT added.
		
		Returns: the ext_node that should be used
		"""
		
		ext_node = "INT" + str(self.internal_nodes)
		self.internal_nodes = self.internal_nodes + 1
		return ext_node
	
	def is_int_node_internal_only(self, int_node):
		"""Returns: 
		True if the node was automatically added by the simulator,
		False, otherwise.
		
		"""
		return self.nodes_dict[int_node].find("INT") > -1
	
	def is_nonlinear(self):
		"""Returns True if at least a element in the circuit is NL.
		"""
		for elem in self.elements:
			if elem.is_nonlinear:
				return True
		return False
	
	def get_locked_nodes(self):
		"""Restituisce una lista che contiene tutti i nodi connessi a elementi non lineari.
		Normalmente questa lista viene poi passata a dc_solve o mdn_solver che la passeranno 
		a get_td.
		"""
		locked_nodes = []
		nl_elements = [elem for elem in self.elements if elem.is_nonlinear]
		#nl_elements = filter(lambda elem: elem.is_nonlinear, element_list)
		for elem in nl_elements:
			oports = elem.get_output_ports()
			for index in range(len(oports)):
				ports = elem.get_drive_ports(index)
				for port in ports:
					locked_nodes.append(port)
		return locked_nodes
	
	def ext_node_to_int(self, ext_node):
		"""This function returns the integer id associated with the external node id, the
		string ext_node.
		
		Parameters:
		ext_node: the external node id to be converted. This is always a string.
		
		Note: this method is slow, because it has to look through circuit_inst.nodes_dict
		
		Throws a NodeNotFoundError exception.

		
		Returns: the int id, int_node associated.
		"""
		items = self.nodes_dict.items()
		values = [value for key, value in items]
		
		try:
			index = values.index(ext_node)
		except ValueError:
			raise NodeNotFoundError, ""
		
		return items[index][0]
		
	def int_node_to_ext(self, int_node):
		"""This function returns the string id associated with the integer internal node id
		int_node.
		
		Parameters:
		int_node: the internal node id to be converted. This is always a integer.
		
		Note: this is the same as circuit_inst.nodes_dict[int_node], except that it throws a
		NodeNotFoundError exception and not a KeyError.
		This method is fast.
		
		Returns: the string id, ext_node associated.
		"""
		try:
			ret = self.nodes_dict[int_node]
		except KeyError:
			raise NodeNotFoundError, ""
		
		return ret

	def has_duplicate_elem(self):
		for index1 in range(len(self.elements)):
			for index2 in range(index1+1, len(self.elements)):
				if self.elements[index1].letter_id == self.elements[index2].letter_id and \
				self.elements[index1].descr == self.elements[index2].descr:
					return True
		return False


	def get_ground_node(self):
		"Returns the (external) reference node. AKA GND."
		return '0'
	def get_elem_by_name(self, name):
		for e in self.elements:
			if e.name == name:
				return e
		return None

	def add_model(self, model_type, model_label, model_parameters):
		"""Add a model to the available models
		Inputs:
		* models (a dictionary, "label":model instance), the available models. None if no model is available/defined.
		* model_type (string), the model type (eg "ekv")
		* model_label (string), a unique identifier for the model being added
		* model_parameters (dict), a dictionary holding the parameters to be 
		supplied to the model to initialize it.

		returns: the updated models
		"""

		if model_type == "ekv":
			model_iter = ekv.ekv_mos_model(**model_parameters)
			model_iter.name = model_label
		elif model_type == "mosq":
			model_iter = mosq.mosq_mos_model(**model_parameters)
			model_iter.name = model_label
		else:
			raise CircuitError, "Unknown model %s" % (model_type,)
		self.models.update({model_label:model_iter})
		return self.models
			
	def remove_model(self, model_label):
		"""Remove a model to the available models
		Inputs:
		* models (a dictionary, "label":model instance), the available models or None if no model is available/defined.
		* model_label (string), the unique identifier corresponding to the model 
		being removed

		This method currently silently ignores models that are not defined.

		returns: None
		"""
		if self.models is not None and self.models.has_key(model_label):
			del self.models[model_label]
		# should print a warning here

	def add_resistor(self, name, ext_n1, ext_n2, R):
		"""Adds a resistor to the circuit (also takes care that the nodes are
		added as well).
	
		Parameters:
		name (string): the resistor name (eg "R1"). The first letter is replaced by an R
		ext_n1, ext_n2 (string): the nodes to which the resistor is connected. 
					eg. "in" or "out_a"
		R (float): resistance (float)
	
		Returns: True
		"""
		n1 = self.add_node_to_circ(ext_n1)
		n2 = self.add_node_to_circ(ext_n2)
	
		if R == 0:
			raise CircuitError, "ZERO-valued resistors are not allowed."

		elem = devices.resistor(n1=n1, n2=n2, R=R)
		elem.descr = name[1:]
		self.elements = self.elements + [elem]
		return True
	
	def add_capacitor(self, name, ext_n1, ext_n2, C, ic=None):
		"""Adds a capacitor to the circuit (also takes care that the nodes are
		added as well).
	
		Parameters:
		name (string): the capacitor name (eg "C1"). The first letter is always C.
		ext_n1, ext_n2 (string): the nodes to which the element is connected. 
					eg. "in" or "out_a"
		C (float): capacitance (float)
		ic (float): initial condition, see simulation types for how this affects
			the results.	

		Returns: True
		"""
		if C == 0:
			raise CircuitError, "ZERO-valued capacitors are not allowed."

		n1 = self.add_node_to_circ(ext_n1)
		n2 = self.add_node_to_circ(ext_n2)
	
		elem = devices.capacitor(n1=n1, n2=n2, C=C, ic=ic)
		elem.descr = name[1:]
	
		self.elements = self.elements + [elem]
		return True

	def add_inductor(self, name, ext_n1, ext_n2, L, ic=None):
		"""Adds an inductor to the circuit (also takes care that the nodes are
		added as well).
	
		Parameters:
		name (string): the inductor name (eg "Lfilter"). The first letter is always L.
		ext_n1, ext_n2 (string): the nodes to which the element is connected. 
					eg. "in" or "out_a"
		C (float): inductance
		ic (float): initial condition, see simulation types for how this affects
			the results.	

		Returns: True
		"""

		n1 = self.add_node_to_circ(ext_n1)
		n2 = self.add_node_to_circ(ext_n2)
	
		elem = devices.inductor(n1=n1, n2=n2, L=L, ic=ic)
		elem.descr = name[1:]
	
		self.elements = self.elements + [elem]
		return True

	def add_inductor_coupling(self, name, L1, L2, Kvalue):
		""" Write DOC XXX
		"""
		L1descr = L1[1:]
		L2descr = L2[1:]
		L1elem, L2elem = None, None

		for e in self.elements:
			if isinstance(e, devices.inductor) and L1descr == e.descr:
				L1elem = e
			elif isinstance(e, devices.inductor) and L2descr == e.descr:
				L2elem = e

		if L1elem is None or L2elem is None:
			error_msg = "One or more coupled inductors for %s were not found: %s (found: %s), %s (found: %s)." % \
			(name, L1, L1elem is not None, L2, L2elem is not None)
			printing.print_general_error(error_msg)
			printing.print_general_error("Quitting.")
			sys.exit(30)

		M = math.sqrt(L1elem.L * L2elem.L) * Kvalue

		elem = devices.inductor_coupling(L1=L1, L2=L2, K=Kvalue, M=M)
		elem.descr = name[1:]
		L1elem.coupling_devices.append(elem)
		L2elem.coupling_devices.append(elem)

		self.elements = self.elements + [elem]

	def add_vsource(self, name, ext_n1, ext_n2, vdc, vac=0, function=None):
		"""Adds a voltage source to the circuit (also takes care that the nodes 
		are added as well).
	
		Parameters:
		name (string): the volatge source name (eg "VA"). The first letter is always V.
		ext_n1, ext_n2 (string): the nodes to which the element is connected. 
					eg. "in" or "out_a"
		vdc (float): DC voltage
		vac (float): AC voltage (optional)
		function (function): optional time function. See devices.py for built-ins.

		Returns: True
		"""
		n1 = self.add_node_to_circ(ext_n1)
		n2 = self.add_node_to_circ(ext_n2)
	
		elem = devices.vsource(n1=n1, n2=n2, vdc=vdc, abs_ac=vac)
		elem.descr = name[1:]
	
		if function is not None:
			elem.is_timedependent = True
			elem._time_function = function
	
		self.elements = self.elements + [elem]
		return True
	
	def add_isource(self, name, ext_n1, ext_n2, idc, iac=0, function=None):
		"""Adds a current source to the circuit (also takes care that the nodes 
		are added as well).
	
		Parameters:
		name (string): the current source name (eg "IA"). The first letter is always I.
		ext_n1, ext_n2 (string): the nodes to which the element is connected. 
					eg. "in" or "out_a"
		idc (float): DC current
		iac (float): AC current (optional)
		function (function): optional time function. See devices.py for built-ins.

		Returns: True
		"""
		n1 = self.add_node_to_circ(ext_n1)
		n2 = self.add_node_to_circ(ext_n2)
	
		elem = devices.isource(n1=n1, n2=n2, idc=idc, abs_ac=iac)
		elem.descr = name[1:]
	
		if function is not None:
			elem.is_timedependent = True
			elem._time_function = function
	
		self.elements = self.elements + [elem]
		return True

	def add_diode(self, name, ext_n1, ext_n2, Is=None, Rs=None, m=None, T=None, ic=None, off=False):
		"""Adds a diode to the circuit (also takes care that the nodes 
		are added as well).
	
		Parameters:
		name (string): the diode name (eg "D1"). The first letter is always D.
		ext_n1, ext_n2 (string): the nodes to which the element is connected. 
					eg. "in" or "out_a"
		Is (float): Inverse sat current (optional)
		rs (float): series resistance (optional)
		m (int): shunt multiplier
		T (float): temperature
		ic (float): initial condition (not implemented yet)

		Returns: True
		"""
		n1 = self.add_node_to_circ(ext_n1)
		n2 = self.add_node_to_circ(ext_n2)
	
		if Rs is not None: #we add a Rs on the anode
			new_node = self.generate_internal_only_node_label()
			self.add_resistor(name="RS_"+name, ext_n1=ext_n1, ext_n2=new_node, R=Rs)
			n1 = self.add_node_to_circ(new_node)
	
		elem = devices.diode(n1=n1, n2=n2, Io=Is, m=m, T=T, ic=ic, off=off)
		elem.descr = name[1:]
		self.elements = self.elements + [elem]

		return True
	
	def add_mos(self, name, ext_nd, ext_ng, ext_ns, ext_nb, w, l, model_label, models, m=None, n=None):
		"""Adds a mosfet to the circuit (also takes care that the nodes 
		are added as well).
	
		Parameters:

		name (string): the mos name (eg "M1"). The first letter is always M.
		ext_nd (string): drain node
		ext_ng (string): gate node
		ext_ns (string): source node
		ext_nb (string): bulk node
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

		nd = self.add_node_to_circ(ext_nd)
		ng = self.add_node_to_circ(ext_ng)
		ns = self.add_node_to_circ(ext_ns)
		nb = self.add_node_to_circ(ext_nb)	

		if not models.has_key(model_label):
			raise ModelError, "Unknown model id: "+model_label
		if isinstance(models[model_label], ekv.ekv_mos_model):
			elem = ekv.ekv_device(nd, ng, ns, nb, w, l, models[model_label], m, n)
		elif isinstance(models[model_label], mosq.mosq_mos_model):
			elem =  mosq.mosq_device(nd, ng, ns, nb, w, l, models[model_label], m, n)
		else:
			raise Exception, "Unknown model type for "+model_label

		#elem = mosq.mosq(nd, ng, ns, kp=kp, w=w, l=l, vt=vt, mos_type=mos_type, lambd=lambd)
		elem.descr = name[1:]
	
		self.elements = self.elements + [elem]
		return True
		
	def add_vcvs(self, name, ext_n1, ext_n2, ext_sn1, ext_sn2, alpha):
		"""Parses a voltage controlled voltage source (vcvs) from the line 
		supplied, adds its nodes to the circuit instance circ and returns a 
		list holding the vcvs element.
	
		Parameters:
		line: the line, if you have already .split()-ed it, set this to None 
		and supply the elements through line_elements.
		circ: the circuit instance.
		line_elements: will be generated by the function from line.split() 
		if set to None.
	
		Returns: [vcvs_elem]
		"""
	
		n1 = self.add_node_to_circ(ext_n1)
		n2 = self.add_node_to_circ(ext_n2)
		sn1 = self.add_node_to_circ(ext_sn1)
		sn2 = self.add_node_to_circ(ext_sn2)
	
		elem = devices.evsource(n1=n1, n2=n2, sn1=sn1, sn2=sn2, alpha=alpha)
		elem.descr = name[1:]
	
		self.elements = self.elements + [elem]
		return True
	
	def add_vccs(self, name, ext_n1, ext_n2, ext_sn1, ext_sn2, alpha):
		"""
		"""

		n1 = self.add_node_to_circ(ext_n1)
		n2 = self.add_node_to_circ(ext_n2)
		sn1 = self.add_node_to_circ(ext_sn1)
		sn2 = self.add_node_to_circ(ext_sn2)
	
		elem = devices.gisource(n1=n1, n2=n2, sn1=sn1, sn2=sn2, alpha=alpha)
		elem.descr = name[1:]
	
		self.elements = self.elements + [elem]
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
			circuit.user_defined_modules_dict.update({module_name:module})

		elem_class = getattr(module, label)
		
		param_dict.update({"convert_units": convert_units})
		param_dict.update({"circuit_node": self.add_node_to_circ})

	
		elem = elem_class(**param_dict)
		elem.descr = name[1:]
		elem.letter_id = "y"
	
		if hasattr(elem, "check"):
			selfcheck_result, error_msg = elem.check()
			if not selfcheck_result:
				raise NetlistParseError, "module: " + module_name + " elem type: "+ elem_type_name+" error: "+\
				error_msg
	
		self.elements = self.elements + [elem]
		return True

	def remove_elem(self, elem):
		"""Removes an element from ther circuit and takes care that no
		"orphan" nodes are left.
		circ: the circuit instance
		elem: the element to be removed

		Returns: True if the element was found and removed, False otherwise
		"""
		if not elem in self.elements:
			return False

		self.elements.remove(elem)
		nodes = []
		if hasattr(elem, n1) and elem.n1 != 0:
			nodes = nodes + [n1]
		if hasattr(elem, n2) and elem.n2 != 0 and not elem.n2 in nodes:
			nodes = nodes + [n2]
		if elem.is_nonlinear:
			for n1, n2 in elem.ports:
				if n1 != 0 and not n1 in nodes:
					nodes = nodes + [n1]
				if n2 != 0 and not n2 in nodes:
					nodes = nodes + [n2]

		remove_list = copy.copy(nodes)
		for n in nodes:
			for e in self.elements:
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

	def find_vde_index(self, id_wdescr):
		vde_index = 0
		for elem in self.elements:
			if is_elem_voltage_defined(elem):
				if (elem.letter_id + elem.descr).upper() == id_wdescr.upper():
					break
				else:
					vde_index += 1
		print "%s found at index %d" % (id_wdescr, vde_index)
		return vde_index

	def find_vde(self, index):
		index = index - len(circ.nodes) + 1
		ni = 0
		found = False
		for e in self.elements:
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
	True se the elem is a vsource, inductor, evsource or hvsource
	False otherwise.
	"""
	if isinstance(elem, devices.vsource) or isinstance(elem, devices.evsource) or \
	isinstance(elem, devices.hvsource) or isinstance(elem, devices.inductor) \
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

class subckt:
	"""This class holds the necessary information about a circuit.
	An instance of this class is returned by 
	  
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
	_before_ calling circ.add_node_to_circ()
	
	While instatiating/calling a subcircuit wrap circ in this.
	"""

	def __init__(self, circ, connection_nodes_dict, subckt_name, subckt_label):
		self.circ = circ
		self.prefix = subckt_label + "-" + subckt_name + "-"
		self.subckt_node_filter_dict = {}
		self.subckt_node_filter_dict.update(connection_nodes_dict)		
		self.subckt_node_filter_dict.update({'0':'0'})
	def add_node_to_circ(self, ext_name):
		"""We want to perform the following conversions:
		connected node in the subcircuit -> node in the upper circuit
		local-only node of the subcircuit -> rename it to something uniq
		REF (0) -> REF (0)
		
		And then call circ.add_node_to_circ()
		"""
		if not self.subckt_node_filter_dict.has_key(ext_name):
			self.subckt_node_filter_dict.update({ext_name:self.prefix+ext_name})
		int_node = self.circ.add_node_to_circ(self.subckt_node_filter_dict[ext_name])
		return int_node
