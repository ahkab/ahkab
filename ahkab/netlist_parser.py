# -*- coding: iso-8859-1 -*-
# netlist_parser.py
# Netlist parser module
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

"""Parse spice-like files, generate circuits instances and print them.
The syntax is explained in the docs and it's based on [1] whenever possible.

Ref. [1] http://newton.ex.ac.uk/teaching/CDHW/Electronics2/userguide/
"""

import sys, imp, math, copy
import ahkab
import circuit
import dc_analysis
import devices, diode, mosq, ekv, switch
import printing, utilities, plotting, options

# analyses syntax
from dc_analysis import specs as dc_spec
from ac import specs as ac_spec
from transient import specs as tran_spec
from pss import specs as pss_spec
from symbolic import specs as symbolic_spec
specs = {}
for i in dc_spec, ac_spec, tran_spec, pss_spec, symbolic_spec:
	specs.update(i)

def parse_circuit(filename, read_netlist_from_stdin=False):
	"""Parse a SPICE-like netlist and return a circuit instance 
	that includes all components, all nodes known
	with that you can recreate mna and N at any time.
	Note that solving the circuit requires accessing to the elements in 
	the circuit instance to evaluate non linear elements' currents.
	
	Directives are collected in a list and returned too, except for
	subcircuits, those are added to circuit.subckts_dict.
	
	Returns:
	(circuit_instance, directives)
	"""
	# Lots of differences with spice's syntax:
	# Support for alphanumeric node names, but the ref has to be 0. always
	# .end is not required, but if is used anything following it is ignored
	# many others, see doc.
	
	circ = circuit.circuit(title="", filename=filename)
	
	if not read_netlist_from_stdin:
		ffile = open(filename, "r")
	else:
		ffile = sys.stdin
		
	file_list = [(ffile, "unknown", not read_netlist_from_stdin)]
	file_index = 0
	directives = []
	model_directives = []
	postproc = []
	subckts_list_temp = []
	netlist_lines = []
	current_subckt_temp = []
	within_subckt = False
	line_n = 0
	
	try:
		while ffile is not None:
			while True:
				line = ffile.readline()
				if len(line) == 0:
					break # check for EOF
				line_n = line_n + 1
				line = line.strip().lower()
				if line_n == 1:
					# the first line is always the title
					circ.title = line
					continue
				elif len(line) == 0:
					continue # empty line is really empty after strip()
				line = join_lines(ffile, line)
				if line[0] == "*": # comments start with *
					continue
				
				# directives are grouped together and evaluated after
				# we have the whole circuit.
				# subcircuits are grouped too, but processed first
				if line[0] == ".":
					line_elements = line.split()
					if line_elements[0] == '.subckt':
						if within_subckt:
							raise NetlistParseError, "nested subcircuit declaration detected"
						current_subckt_temp = current_subckt_temp + [(line, line_n)]
						within_subckt = True
					elif line_elements[0] == '.ends':
						if not within_subckt:
							raise NetlistParseError, "corresponding .subckt not found"
						within_subckt = False
						subckts_list_temp.append(current_subckt_temp)
						current_subckt_temp = []
					elif line_elements[0] == '.include':
						file_list.append(parse_include_directive(line, line_elements=None))
					elif line_elements[0] == ".end":
						break
					elif line_elements[0] == ".plot":
						postproc.append((line, line_n))
					elif line_elements[0] == ".model":
						model_directives.append((line, line_n))
					else:
						directives.append((line, line_n))
					continue
				
				if within_subckt:
					current_subckt_temp  = current_subckt_temp + [(line, line_n)]
				else:
					netlist_lines = netlist_lines + [(line, line_n)]
			if within_subckt:
				raise NetlistParseError, ".ends not found"
			
			file_index = file_index + 1
			ffile = get_next_file_and_close_current(file_list, file_index)
			#print file_list
			
	except NetlistParseError, (msg,):
		if len(msg):
			printing.print_general_error(msg)
		printing.print_parse_error(line_n, line)
		#if not read_netlist_from_stdin:
			#ffile.close()
		sys.exit(45)
	
	#if not read_netlist_from_stdin:
		#ffile.close()

	models = parse_models(model_directives)
	
	# now we parse the subcircuits, we want a circuit.subckt object that holds the netlist code,
	# the nodes and the subckt name in a handy way.
	# We will create all the elements in the subckt every time it is instantiated in the netlist file.
	subckts_dict = {}
	for subckt_temp in subckts_list_temp:
		subckt_obj = parse_sub_declaration(subckt_temp)
		if not subckts_dict.has_key(subckt_obj.name):
			subckts_dict.update({subckt_obj.name:subckt_obj})
		else:
			raise NetlistParseError, "subckt " + subckt_obj.name + " has been redefined"

	elements = main_netlist_parser(circ, netlist_lines, subckts_dict, models)
	circ.elements = circ.elements + elements
	
	return (circ, directives, postproc)
	
def main_netlist_parser(circ, netlist_lines, subckts_dict, models):
	elements = []
	try:
		for (line, line_n) in netlist_lines:
			# elements: detect the element type and call the
			# appropriate parsing function
			# we always use normal convention V opposite to I
			# n1 is +, n2 is -, current flows from + to -
			line_elements = line.split()
			if line[0] == "r":
				elements = elements + \
				parse_elem_resistor(line, circ, line_elements)
			elif line[0] == "c":
				elements = elements + \
				parse_elem_capacitor(line, circ, line_elements)
			elif line[0] == "l":
				elements = elements + \
				parse_elem_inductor(line, circ, line_elements)
			elif line[0] == "k":
				elements = elements + \
				parse_elem_inductor_coupling(line, circ, \
				line_elements, elements)
			elif line[0] == "v":
				elements = elements + \
				parse_elem_vsource(line, circ, line_elements)
			elif line[0] == "i":
				elements = elements + \
				parse_elem_isource(line, circ, line_elements)
			elif line[0] == "d":
				elements = elements + \
				parse_elem_diode(line, circ, line_elements, models)
			elif line[0] == 'm': #mosfet
				elements = elements + \
				parse_elem_mos(line, circ, line_elements, models)
			elif line[0] == "e": #vcvs
				elements = elements + \
				parse_elem_vcvs(line, circ, line_elements)
			elif line[0] == "g": #vccs
				elements = elements + \
				parse_elem_vccs(line, circ, line_elements)
			elif line[0] == "s": #vccs
				elements = elements + \
				parse_elem_switch(line, circ, line_elements, models)
			elif line[0] == "y": #User defined module -> MODIFY
				elements = elements + \
				parse_elem_user_defined(line, circ, line_elements)
			elif line[0] == "x": #User defined module -> MODIFY
				elements = elements + \
				parse_sub_instance(line, circ, subckts_dict, line_elements, models)
			else:
				raise NetlistParseError, "unknown element."
	except NetlistParseError, (msg,):
		if len(msg):
			printing.print_general_error(msg)
		printing.print_parse_error(line_n, line)
		sys.exit(45)
	
	return elements

def get_next_file_and_close_current(file_list, file_index):
	if file_list[file_index - 1][2]:
		file_list[file_index - 1][0].close()
	if file_index == len(file_list):
		ffile = None
	else:
		ffile = open(file_list[file_index][1], "r")
		file_list[file_index][0] = ffile
	return ffile

def parse_models(models_lines):
	models = {}
	for line, line_n in models_lines:
		tokens = line.replace("(","").replace(")", "").split()
		if len(tokens) < 3:
			raise NetlistParseError, ("Syntax error in model declaration on line " + str(line_n) + ".\n\t"+line,)
		model_label = tokens[2]
		model_type = tokens[1]
		model_parameters = {}
		for index in range(3, len(tokens)):
			if tokens[index][0] == "*":
				break
			(label, value) = parse_param_value_from_string(tokens[index])
			model_parameters.update({label.upper():value})
		if model_type == "ekv":
			model_iter = ekv.ekv_mos_model(**model_parameters)
			model_iter.name = model_label
		elif model_type == "mosq":
			model_iter = mosq.mosq_mos_model(**model_parameters)
			model_iter.name = model_label
		elif model_type == "diode" or model_type == 'd':
			model_parameters.update({'name':model_label})
			model_iter = diode.diode_model(**model_parameters)
		elif model_type == "sw":
			model_parameters.update({'name':model_label})
			model_iter = switch.vswitch_model(**model_parameters)
		#elif model_type == "csw":
		#	model_parameters.update({'name':model_label})
		#	model_iter = switch.iswitch_model(**model_parameters)
		else:
			raise NetlistParseError, ("Unknown model ("+model_type+") on line " + str(line_n) + ".\n\t"+line,)
		models.update({model_label:model_iter})
	return models
			

def parse_elem_resistor(line, circ, line_elements=None):
	"""Parses a resistor from the line supplied, adds its nodes to the circuit
	instance circ and returns a list holding the resistor element.
	
	Parameters:
	line: the line, if you have already .split()-ed it, set this to None 
	and supply the elements through line_elements.
	circ: the circuit instance.
	line_elements: will be generated by the function from line.split() 
	if set to None.
	
	Returns: [resistor_elem]
	"""
	if line_elements is None:
		line_elements = line.split()
	
	if (len(line_elements) < 4) or (len(line_elements) > 4 and not line_elements[4][0] == "*"):
		raise NetlistParseError, ""

	ext_n1 = line_elements[1]
	ext_n2 = line_elements[2]
	n1 = circ.add_node(ext_n1)
	n2 = circ.add_node(ext_n2)
	
	R = convert_units(line_elements[3])

	if R == 0:
		raise NetlistParseError, "ZERO-valued resistors are not allowed." 

	elem = devices.resistor(n1=n1, n2=n2, R=R)
	elem.descr = line_elements[0][1:]
	
	return [elem]
	
def parse_elem_capacitor(line, circ, line_elements=None):
	"""Parses a capacitor from the line supplied, adds its nodes to the circuit
	instance circ and returns a list holding the capacitor element.
	
	Parameters:
	line: the line, if you have already .split()-ed it, set this to None 
	and supply the elements through line_elements.
	circ: the circuit instance.
	line_elements: will be generated by the function from line.split() 
	if set to None.
	
	Returns: [capacitor_elem]
	"""
	if line_elements is None:
		line_elements = line.split()
		
	if (len(line_elements) < 4) or (len(line_elements) > 5 and not line_elements[6][0]=="*"):
		raise NetlistParseError, ""
		
	ic = None
	if len(line_elements) == 5 and not line_elements[4][0] == '*':
		(label, value) = parse_param_value_from_string(line_elements[4])
		if label == "ic":
			ic = convert_units(value)
		else:
			raise NetlistParseError, "unknown parameter " + label

	ext_n1 = line_elements[1]
	ext_n2 = line_elements[2]
	n1 = circ.add_node(ext_n1)
	n2 = circ.add_node(ext_n2)
	
	elem = devices.capacitor(n1=n1, n2=n2, C=convert_units(line_elements[3]), ic=ic)
	elem.descr = line_elements[0][1:]
	
	return [elem]

def parse_elem_inductor(line, circ, line_elements=None):
	"""Parses a inductor from the line supplied, adds its nodes to the circuit
	instance circ and returns a list holding the inductor element.
	
	Parameters:
	line: the line, if you have already .split()-ed it, set this to None 
	and supply the elements through line_elements.
	circ: the circuit instance.
	line_elements: will be generated by the function from line.split() 
	if set to None.
	
	Returns: [inductor_elem]
	"""
	if line_elements is None:
		line_elements = line.split()
	
	if (len(line_elements) < 4) or (len(line_elements) > 5 and not line_elements[6][0]=="*"):
		raise NetlistParseError, ""
		
	ic = None
	if len(line_elements) == 5 and not line_elements[4][0] == '*':
		(label, value) = parse_param_value_from_string(line_elements[4])
		if label == "ic":
			ic = convert_units(value)
		else:
			raise NetlistParseError, "unknown parameter " + label
		
	ext_n1 = line_elements[1]
	ext_n2 = line_elements[2]
	n1 = circ.add_node(ext_n1)
	n2 = circ.add_node(ext_n2)
	
	elem = devices.inductor(n1=n1, n2=n2, L=convert_units(line_elements[3]), ic=ic)
	elem.descr = line_elements[0][1:]
	
	return [elem]

def parse_elem_inductor_coupling(line, circ, line_elements=None, elements=[]):
	"""Parses a inductor coupling from the line supplied, 
	returns a list holding the element.
	
	Parameters:
	line: the line, if you have already .split()-ed it, set this to None 
	and supply the elements through line_elements.
	circ: the circuit instance.
	line_elements: will be generated by the function from line.split() 
	if set to None.
	
	Returns: [inductor_coupling_elem]
	"""
	if line_elements is None:
		line_elements = line.split()
	
	if (len(line_elements) < 4) or (len(line_elements) > 4 and not line_elements[5][0]=="*"):
		raise NetlistParseError, ""
	
	name = line_elements[0]
	L1 = line_elements[1]
	L2 = line_elements[2]

	try:
		Kvalue = convert_units(line_elements[3])
	except ValueError:
		(label, value) = parse_param_value_from_string(line_elements[3])
		if not label == "k":
			raise NetlistParseError, "unknown parameter " + label
		Kvalue = convert_units(value)

	L1descr = L1[1:]
	L2descr = L2[1:]
	L1elem, L2elem = None, None

	for e in elements:
		if isinstance(e, devices.inductor) and L1descr == e.descr:
			L1elem = e
		elif isinstance(e, devices.inductor) and L2descr == e.descr:
			L2elem = e

	if L1elem is None or L2elem is None:
		error_msg = "One or more coupled inductors for %s were not found: %s (found: %s), %s (found: %s)." % \
		(name, L1, L1elem is not None, L2, L2elem is not None)
		printing.print_general_error(error_msg)
		raise NetlistParseError, ""

	M = math.sqrt(L1elem.L * L2elem.L) * Kvalue

	elem = devices.inductor_coupling(L1=L1, L2=L2, K=Kvalue, M=M)
	elem.descr = name[1:]
	L1elem.coupling_devices.append(elem)
	L2elem.coupling_devices.append(elem)	
	
	return [elem]

def parse_elem_vsource(line, circ, line_elements=None):
	"""Parses a vsource from the line supplied, adds its nodes to the circuit
	instance circ and returns a list holding the vsource element.
	
	Parameters:
	line: the line, if you have already .split()-ed it, set this to None 
	and supply the elements through line_elements.
	circ: the circuit instance.
	line_elements: will be generated by the function from line.split() 
	if set to None.
	
	Returns: [vsource_elem]
	"""
	if line_elements is None:
		line_elements = line.split()
	
	if len(line_elements) < 3:
		raise NetlistParseError, ""
	
	vdc = None
	vac = None
	function = None
	
	index = 3
	while True: #for index in range(3, len(line_elements)):
		if index == len(line_elements):
			break
		if line_elements[index][0] == '*':
			break
		
		(label, value) = parse_param_value_from_string(line_elements[index])
		
		if label == 'type':
			if value == 'vdc':
				param_number = 0
			elif value == 'vac':
				param_number = 0
			elif value == 'pulse':
				param_number = 7
			elif value == 'exp':
				param_number = 6
			elif value == 'sin':
				param_number = 5
			else:
				raise NetlistParseError, "unknown signal type."
			if param_number and function is None:
				function = parse_time_function(value, line_elements[index+1:index+param_number+1], "voltage")
				index = index + param_number
				#continue
			elif function is not None:
				raise NetlistParseError, "only a time function can be defined."
		elif label == 'vdc':
			vdc = convert_units(value)
		elif label == 'vac':
			vac = convert_units(value)
		else:
			raise NetlistParseError, ""
		index = index + 1

	if vdc == None and function == None:
		raise NetlistParseError, "neither vdc nor a time function are defined."
	
	#usual
	ext_n1 = line_elements[1]
	ext_n2 = line_elements[2]
	n1 = circ.add_node(ext_n1)
	n2 = circ.add_node(ext_n2)
	
	elem = devices.vsource(n1=n1, n2=n2, vdc=vdc, abs_ac=vac)
	elem.descr = line_elements[0][1:]
	
	if function is not None:
		elem.is_timedependent = True
		elem._time_function = function
	
	return [elem]
	
def parse_elem_isource(line, circ, line_elements=None):
	"""Parses a isource from the line supplied, adds its nodes to the circuit
	instance circ and returns a list holding the isource element.
	
	Parameters:
	line: the line, if you have already .split()-ed it, set this to None 
	and supply the elements through line_elements.
	circ: the circuit instance.
	line_elements: will be generated by the function from line.split() 
	if set to None.
	
	Returns: [isource_elem]
	"""
	if line_elements is None:
		line_elements = line.split()
	
	if len(line_elements) < 3:
		raise NetlistParseError, ""
	
	idc = None
	iac = None
	function = None
	
	index = 3
	while True: #for index in range(3, len(line_elements)):
		if index == len(line_elements):
			break
		if line_elements[index][0] == '*':
			break
		
		(label, value) = parse_param_value_from_string(line_elements[index])
		
		if label == 'type':
			if value == 'idc':
				param_number = 0
			elif value == 'iac':
				param_number = 0
			elif value == 'pulse':
				param_number = 7
			elif value == 'exp':
				param_number = 6
			elif value == 'sin':
				param_number = 5
			else:
				raise NetlistParseError, "unknown signal type."
			if param_number and function is None:
				function = parse_time_function(value, line_elements[index+1:index+param_number+1], "current")
				index = index + param_number
			elif function is not None:
				raise NetlistParseError, "only a time function can be defined."
		elif label == 'idc':
			idc = convert_units(value)
		elif label == 'iac':
			iac = convert_units(value)
		else:
			raise NetlistParseError, ""
		index = index + 1

	if idc == None and function == None:
		raise NetlistParseError, "neither idc nor a time function are defined."
	
	ext_n1 = line_elements[1]
	ext_n2 = line_elements[2]
	n1 = circ.add_node(ext_n1)
	n2 = circ.add_node(ext_n2)
	
	elem = devices.isource(n1=n1, n2=n2, idc=idc, abs_ac=iac)
	elem.descr = line_elements[0][1:]
	
	if function is not None:
		elem.is_timedependent = True
		elem._time_function = function
	
	return [elem]

def parse_elem_diode(line, circ, line_elements=None, models=None):
	"""Parses a diode from the line supplied, adds its nodes to the circuit
	instance circ and returns a list holding the diode element and a resistor,
	if the diode has Rs != 0.
	
	Diode syntax:
	#DX N+ N- <MODEL_LABEL> <AREA=xxx>
	
	Parameters:
	line: the line, if you have already .split()-ed it, set this to None 
	and supply the elements through line_elements.
	circ: the circuit instance.
	line_elements: will be generated by the function from line.split() 
	if set to None.
	
	Returns: a list with the diode element and a optional resistor Rs
	"""
	#sarebbe bello implementare anche: <IC=VD> <TEMP=T>
	if line_elements is None:
		line_elements = line.split()
	
	Area = None
	T  = None
	ic = None
	off = False
	
	if (len(line_elements) < 4):
		raise NetlistParseError, ""
	
	model_label = line_elements[3] 

	for index in range(4, len(line_elements)):
		if line_elements[index][0] == '*':
			break
		(param, value) = parse_param_value_from_string(line_elements[index])
		
		value = convert_units(value)
		if param == "area":
			Area = value
		elif param == "t":
			T = value
		elif param == "ic":
			ic = value
		elif param == "off":
			if not len(value):
				off = True
			else:
				off = convert_boolean(value)
		else:
			raise NetlistParseError, "unknown parameter " + param
		
	ext_n1 = line_elements[1]
	ext_n2 = line_elements[2]
	n1 = circ.add_node(ext_n1)
	n2 = circ.add_node(ext_n2)
	
	#if Rs: #we need to add a Rs on the anode
	#	new_node = n1
	#	n1 = circ.generate_internal_only_node_label()
	#	#print "-<<<<<<<<"+str(n1)+" "+str(n2) +" "+str(new_node)
	#	rs_elem = devices.resistor(n1=new_node, n2=n1, R=Rs)
	#	rs_elem.descr = "INT"
	#	return_list = return_list + [rs_elem]
	
	if not models.has_key(model_label):
		raise NetlistParseError, "Unknown model id: "+model_label
	elem = diode.diode(n1=n1, n2=n2, model=models[model_label], AREA=Area, T=T, ic=ic, off=off)
	elem.descr = line_elements[0][1:]
	return [elem]

def parse_elem_mos(line, circ, line_elements, models):
	"""Parses a mos from the line supplied, adds its nodes to the circuit
	instance circ and returns a list holding the mos element.
	
	MOS syntax:
	MX ND NG NS KP=xxx Vt=xxx W=xxx L=xxx type=n/p <LAMBDA=xxx>
	
	Parameters:
	line: the line, if you have already .split()-ed it, set this to None 
	and supply the elements through line_elements.
	circ: the circuit instance.
	line_elements: will be generated by the function from line.split() 
	if set to None.
	
	Returns: [mos_elem]
	"""
	#sarebbe bello implementare anche: <IC=VD> <TEMP=T>
	#(self, nd, ng, ns, kp=1.0e-14, w, l, mos_type='n', lambd=0, type_of_elem="mosq")
	if line_elements is None:
		line_elements = line.split()

	if (len(line_elements) < 6):
		raise NetlistParseError, "required parameters are missing."
		#print "MX ND NG NS model_id W=xxx L=xxx"
	
	model_label = line_elements[5]
	
	#kp = None
	w  = None
	l = None
	#mos_type = None
	#vt = None
	m = 1
	n = 1		
	#lambd = 0 # va is supposed infinite if not specified
	for index in range(6, len(line_elements)):
		if line_elements[index][0] == '*':
			break
		
		(param, value) = parse_param_value_from_string(line_elements[index])
		
		#if param == 'kp':
		#	kp = convert_units(value)
		if param == "w":
			w = convert_units(value)
		elif param == "l":
			l = convert_units(value)
		#elif param == "vt":
		#	vt = convert_units(value)
		elif param == "m":
			m = convert_units(value)
		#elif param == "vt":
		#	n = convert_units(value)
		#elif param == "type":
		#	if value != 'n' and value != 'p':
		#		raise NetlistParseError, "unknown mos type "+value
		#	mos_type = value
		else:
			raise NetlistParseError, "unknown parameter " + param

	if (w is None) or (l is None):
		raise NetlistParseError, "required parameter is missing."
		#print "MX ND NG NS W=xxx L=xxx <LAMBDA=xxx>"
		
	ext_nd = line_elements[1]
	ext_ng = line_elements[2]
	ext_ns = line_elements[3]
	ext_nb = line_elements[4]
	nd = circ.add_node(ext_nd)
	ng = circ.add_node(ext_ng)
	ns = circ.add_node(ext_ns)
	nb = circ.add_node(ext_nb)	

	if not models.has_key(model_label):
		raise NetlistParseError, "Unknown model id: "+model_label
	if isinstance(models[model_label], ekv.ekv_mos_model):
		elem = ekv.ekv_device(nd, ng, ns, nb, w, l, models[model_label], m, n)
	elif isinstance(models[model_label], mosq.mosq_mos_model):
		elem = mosq.mosq_device(nd, ng, ns, nb, w, l, models[model_label], m, n)
	else:
		raise NetlistParseError, "Unknown MOS model type: "+model_label

	elem.descr = line_elements[0][1:]
	
	return [elem]
		
def parse_elem_vcvs(line, circ, line_elements=None):
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
	if line_elements is None:
		line_elements = line.split()
	
	if (len(line_elements) < 6) or (len(line_elements) > 6 and not line_elements[6][0] == "*"):
		raise NetlistParseError, ""
	
	ext_n1 = line_elements[1]
	ext_n2 = line_elements[2]
	ext_sn1 = line_elements[3]
	ext_sn2 = line_elements[4]
	n1 = circ.add_node(ext_n1)
	n2 = circ.add_node(ext_n2)
	sn1 = circ.add_node(ext_sn1)
	sn2 = circ.add_node(ext_sn2)
	
	elem = devices.evsource(n1=n1, n2=n2, sn1=sn1, sn2=sn2, alpha=convert_units(line_elements[5]))
	elem.descr = line_elements[0][1:]
	
	return [elem]
	
def parse_elem_vccs(line, circ, line_elements=None):
	"""Parses a voltage controlled current source (vccs) from the line 
	supplied, adds its nodes to the circuit instance circ and returns a 
	list holding the vccs element.
	
	Syntax:
	GX N+ N- NC+ NC- VALUE	
	
	Parameters:
	line: the line, if you have already .split()-ed it, set this to None 
	and supply the elements through line_elements.
	circ: the circuit instance.
	line_elements: will be generated by the function from line.split() 
	if set to None.
	
	Returns: [vccs_elem]
	"""

	if line_elements is None:
		line_elements = line.split()
	
	if (len(line_elements) < 6) or (len(line_elements) > 6 \
	and not line_elements[6][0]=="*"):
		raise NetlistParseError, ""
	
	ext_n1 = line_elements[1]
	ext_n2 = line_elements[2]
	ext_sn1 = line_elements[3]
	ext_sn2 = line_elements[4]
	n1 = circ.add_node(ext_n1)
	n2 = circ.add_node(ext_n2)
	sn1 = circ.add_node(ext_sn1)
	sn2 = circ.add_node(ext_sn2)
	
	elem = devices.gisource(n1=n1, n2=n2, sn1=sn1, sn2=sn2, alpha=convert_units(line_elements[5]))
	elem.descr = line_elements[0][1:]
	
	return [elem]

def parse_elem_switch(line, circ, line_elements=None, models=None):
	"""Parses a switch device from the line supplied, adds its nodes to 
	the circuit instance circ and returns a list holding the switch element.
	
	Parameters:
	line: the line, if you have already .split()-ed it, set this to None 
	and supply the elements through line_elements.
	circ: the circuit instance.
	line_elements: will be generated by the function from line.split() 
	if set to None.
	models: the dict of the currently defined models.
	
	General syntax:
	SW1 n1 n2 ns1 ns2 model_label

	Returns: [switch_elem]
	"""
	if line_elements is None:
		line_elements = line.split()
	
	if (len(line_elements) < 6) or (len(line_elements) > 6 and not line_elements[6][0] == "*"):
		raise NetlistParseError, ""
	
	ext_n1 = line_elements[1]
	ext_n2 = line_elements[2]
	ext_sn1 = line_elements[3]
	ext_sn2 = line_elements[4]
	n1 = circ.add_node(ext_n1)
	n2 = circ.add_node(ext_n2)
	sn1 = circ.add_node(ext_sn1)
	sn2 = circ.add_node(ext_sn2)
	
	model_label = line_elements[5]
	if not models.has_key(model_label):
		raise NetlistParseError, "Unknown model id: "+model_label
	if isinstance(models[model_label], switch.vswitch_model):
		elem = switch.switch_device(n1, n2, sn1, sn2, models[model_label])
	#if isinstance(models[model_label], switch.vswitch_model) or \
	#   isinstance(models[model_label], switch.iswitch_model):
	#	elem = mosq.mosq_device(nd, ng, ns, nb, w, l, models[model_label], m, n)
	else:
		raise NetlistParseError, "Unknown MOS model type: "+model_label

	elem.descr = line_elements[0][1:]
	
	return [elem]
	
def parse_elem_user_defined(line, circ, line_elements=None):
	"""Parses a user defined element.
	
	In order for this to work, you should write a module that supplies the
	elem class.
	
	Syntax:
	Y<X> <n1> <n2> module=<module_name> type=<type> [<param1>=<value1> ...]
	
	This method will attempt to load the module <module_name> and it will
	then look for a class named <type>.
	
	An object will be instatiated with the following arguments: 
	n1, n2, param_dict, get_int_id_func, convert_units_func
	Where:
	n1: is the anode of the element
	n2: is the cathode
	param_dict: is a dictionary, its elements are {param1:value1, ...}
	get_int_id_func, convert_units_func are two function that may be used
	in the __init__ method, if needed.
	get_int_id_func: a function that gives back the internal name of a node
	convert_units_func: utility function to convert eg 1p -> 1e-12
	
	See ideal_oscillators.py for a reference implementation.
	
	Parameters:
	line: the line, if you have already .split()-ed it, set this to None 
	and supply the elements through line_elements.
	circ: the circuit instance.
	line_elements: will be generated by the function from line.split() 
	if set to None.
	
	Returns: [userdef_elem]
	"""
	if line_elements is None:
		line_elements = line.split()
	
	if len(line_elements) < 4:
		raise NetlistParseError, ""
	
	param_dict = {}
	for index in range(3, len(line_elements)):
		if line_elements[index][0] == '*':
			break
		
		(param, value) = parse_param_value_from_string(line_elements[index])
		
		if not param_dict.has_key(param):
			param_dict.update({param:value})
		else:
			raise NetlistParseError, param+" already defined."
		
	if param_dict.has_key("module"):
		module_name = param_dict.pop("module", None)
	else:
		raise NetlistParseError, "module name is missing."
	
	if circuit.user_defined_modules_dict.has_key(module_name):
		module = circuit.user_defined_modules_dict[module_name]
	else:
		try:
			fp, pathname, description = imp.find_module(module_name)
			module = imp.load_module(module_name, fp, pathname, description)
		except ImportError:
			raise NetlistParseError, "module " + module_name + " not found."
		circuit.user_defined_modules_dict.update({module_name:module})
	
	if param_dict.has_key("type"):
		elem_type_name = param_dict.pop("type", None)
	else:
		raise NetlistParseError, "type of element is missing."
	
	try:
		elem_class = getattr(module, elem_type_name)
	except AttributeError:
		raise NetlistParseError, "module doesn't have elem type: "+ elem_type_name
		
	ext_n1 = line_elements[1]
	ext_n2 = line_elements[2]
	n1 = circ.add_node(ext_n1)
	n2 = circ.add_node(ext_n2)
	
	elem = elem_class(n1, n2, param_dict, circ.add_node, convert_units)
	elem.descr = line_elements[0][1:]
	elem.letter_id = "y"
	
	selfcheck_result, error_msg = elem.check()
	if not selfcheck_result:
		raise NetlistParseError, "module: " + module_name + " elem type: "+ elem_type_name+" error: "+\
		error_msg
	#fixme non so sicuro che sia una buona idea
	
	return [elem]
	

def parse_time_function(ftype, line_elements, stype):
	"""Parses a time function of type ftype from the line_elements supplied.
	
	ftype: a string, one among "pulse", "exp", "sin"
	line_elements: mustn't hold the "type=<ftype>" element
	stype: set this to "current" for current sources, "voltage" for voltage sources
	
	See devices.pulse, devices.sin, devices.exp for more.
	
	Returns: a time-<function instance
	"""
	if ftype == 'pulse':
		function = parse_pulse_time_function(line_elements, stype)
	elif ftype == 'exp':
		function = parse_exp_time_function(line_elements, stype)
	elif ftype == 'sin':
		function = parse_sin_time_function(line_elements, stype)
	else:
		raise NetlistParseError, "unknown signal type."
	# This way it knows if it is a v/i source in __str__
	function._type = "V"*(stype.lower()=="voltage") + "I"*(stype.lower()=="current")
	return function

def parse_pulse_time_function(line_elements, stype):
	"""This is called by parse_time_function() to actually parse this
	type of functions.
	"""
	function = devices.pulse()
	for token in line_elements:
		(param, value) = parse_param_value_from_string(token)
		if stype == "voltage" and param == 'v1' or stype == "current" \
		and param == 'i1':
			function.v1 = convert_units(value)
		elif stype == "voltage" and param == 'v2' or stype == "current" \
		and param == 'i2':
			function.v2 = convert_units(value)
		elif param == 'td':
			function.td = convert_units(value)
		elif param == 'per':
			function.per = convert_units(value)
		elif param == 'tr':
			function.tr = convert_units(value)
		elif param == 'tf':
			function.tf = convert_units(value)
		elif param == 'pw':
			function.pw = convert_units(value)
		else:
			raise NetlistParseError, "unknown param for time function pulse: "+param
			
	if not function.ready():
		raise NetlistParseError, \
		"required parameters are missing for time function pulse."
	return function
		
def parse_exp_time_function(line_elements, stype):
	"""This is called by parse_time_function() to actually parse this
	type of functions.
	"""
	function = devices.exp()
	for token in line_elements:
		(param, value) = parse_param_value_from_string(token)
		if stype == "voltage" and param == 'v1' or \
		stype == "current" and param == 'i1':
			function.v1 = convert_units(value)
		elif stype == "voltage" and param == 'v2' or \
		stype == "current" and param == 'i2':
			function.v2 = convert_units(value)
		elif param == 'td1':
			function.td1 = convert_units(value)
		elif param == 'tau1':
			function.tau1 = convert_units(value)
		elif param == 'td2':
			function.td2 = convert_units(value)
		elif param == 'tau2':
			function.tau2 = convert_units(value)
		else:
			raise NetlistParseError, "unknown param for time function exp: "+param
			
	if not function.ready():
		raise NetlistParseError, "required params are missing for time function exp."
	return function

def parse_sin_time_function(line_elements, stype):
	"""This is called by parse_time_function() to actually parse this
	type of functions.
	"""
	function = devices.sin()
	for token in line_elements:
		(param, value) = parse_param_value_from_string(token)
		if stype == "voltage" and param == 'vo' \
		or stype == "current" and param == 'io':
			function.vo = convert_units(value)
		elif stype == "voltage" and param == 'va' or \
		stype == "current" and param == 'ia':
			function.va = convert_units(value)
		elif param == 'freq':
			function.freq = convert_units(value)
		elif param == 'theta':
			function.theta = convert_units(value)
		elif param == 'td':
			function.td = convert_units(value)
		else:
			raise NetlistParseError, "unknown param for time function sin: "+param
			
	if not function.ready():
		raise NetlistParseError, "required params are missing for time function sin."
	return function

def convert_units(string_value):
	"""Converts a value conforming to spice's syntax to float. 
	Quote from spice3's manual:
	A number field may be an integer field (eg 12, -44), a floating point 
	field (3.14159), either an integer or floating point number followed by
	an integer exponent (1e-14, 2.65e3), or either an integer or a floating
	point number followed by one of the following scale factors:
	T = 1e12, G = 1e9, Meg = 1e6, K = 1e3, mil = 25.4x1e-6, m = 1e-3,
	u = 1e-6, n = 1e-9, p = 1e-12, f = 1e-15
	
	Raises ValueError if the supplied string can't be interpreted according 
	to the above.
	
	Returns a float.
	"""

	if type(string_value) is float:
		return string_value # not actually a string!
	if not len(string_value):
		raise NetlistParseError, ""

	index = 0
	string_value = string_value.strip().upper()
	while(True):
		if len(string_value)==index:
			break
		if not (string_value[index].isdigit() or string_value[index]=="." or \
			string_value[index]=="+" or  string_value[index]=="-" or \
			string_value[index]=="E"):
			break
		index = index+1
	if index == 0:
		#print string_value
		raise ValueError, "Unable to parse value: %s" % string_value
		#return 0
	numeric_value = float(string_value[:index])
	multiplier = string_value[index:]
	if len(multiplier)==0:
		pass
		#return numeric_value
	elif len(multiplier)==1:
		if multiplier == "T":
			numeric_value = numeric_value*1e12
		elif multiplier == "G":
			numeric_value = numeric_value*1e9
		elif multiplier == "K":
			numeric_value = numeric_value*1e3
		elif multiplier == "M":
			numeric_value = numeric_value*1e-3
		elif multiplier == "U":
			numeric_value = numeric_value*1e-6
		elif multiplier == "N":
			numeric_value = numeric_value*1e-9
		elif multiplier == "P":
			numeric_value = numeric_value*1e-12
		elif multiplier == "F":
			numeric_value = numeric_value*1e-15
		else:
			raise ValueError
	elif len(multiplier)==3:
		if multiplier == "MEG":
			numeric_value = numeric_value*1e6
		elif multiplier == "MIL":
			numeric_value = numeric_value*25.4e-6
		else:
			raise ValueError
	else:
		raise ValueError
	return numeric_value

def parse_postproc(circ, postproc_direc):
	postproc_list = []
	for line, line_n in postproc_direc:
		if line[0] == ".":
			try:
				line_elements = line.split()
				#plot
				if line_elements[0] == ".plot":
					plot_postproc = {}
					plot_postproc["type"] = "plot"
					plot_postproc["analysis"] = line_elements[1]
					if plot_postproc["analysis"] == "tran" or plot_postproc["analysis"] == "shooting":  	
						plot_postproc["x"] = "T"
					elif plot_postproc["analysis"] == "ac":
						plot_postproc["x"] = "w" #QUIRKY FOR THE TIME BEING
					elif plot_postproc["analysis"] == "dc":
						plot_postproc["x"] = an["source_name"]
					else:
						printing.print_general_error("Plotting is unsupported for analysis type "+plot_postproc["analysis"])
						
					graph_labels = ""
					for glabel in line_elements[2:]:
						graph_labels = graph_labels + " " + glabel
					
					l2l1 = plotting.split_netlist_label(graph_labels)

					if plot_postproc["analysis"] == "ac":
						l2l1ac = []
						for l2, l1 in l2l1:
							if l1 is not None:
								l1 = "|%s|" % (l1, )
							else:
								l1 = None
							if l2 is not None:
								l2 = "|%s|" % (l2, )
							else:
								l2 = None
							l2l1ac.append((l2,l1))
						l2l1 = l2l1ac				
					plot_postproc["l2l1"] = l2l1
					postproc_list.append(plot_postproc)
				else:
					raise NetlistParseError("Unknown postproc directive.")
			except NetlistParseError, (msg,):
				if len(msg):
					printing.print_general_error(msg)
				printing.print_parse_error(line_n, line)
				sys.exit(0)
	return postproc_list
					
def parse_ics(directives):
	ics = []
	for line, line_n in directives:
		if line[0] != '.':
			continue
		if line[:3] == '.ic':
			ics += [parse_ic_directive(line)]
	return ics

def parse_analysis(circ, directives):
	"""Parses the analyses.

	Parameters:
	circ: a circuit class instance that descirbes the circuit.
	directives: a list of tuples: (line, line_number). Those lines are taken 
	from the netlist and are the ones that hold the information about the 
	simulations to be performed.
	
	Both of them are returned by parse_circuit()
	
	Returns:
	a list of the analysis, see the code.
	"""
	for line, line_n in directives:
		if line[0] != '.' or line[:3] == '.ic':
			continue
		line_elements = line.split()
		yield parse_single_analysis(line, line_elements)
	return 

def parse_temp_directive(line, line_elements=None):
	"""Parses a TEMP directive:
	
	The syntax is:
	.TEMP <VALUE>
	"""
	if line_elements is None:
		line_elements = line.split()
	
	for token in line_elements[1:]:
		if token[0] == "*":
			break
		value = convert_units(token)
		
	return {"type":"temp", "temp":value}
	
def parse_single_analysis(line, line_elements=None):
	"""Parses an analysis and returns a dict with its parameters.
	"""
	if line_elements is None:
		line_elements = line.split()
	
	an_type = line_elements[0].replace(".", "").lower()
	if not an_type in specs:
		raise NetlistParseError, "Unknown directive: %s" % an_type
	params = list(copy.deepcopy(specs[an_type]['tokens']))
	
	an = {'type':an_type}
	for i in range(len(line_elements[1:])):
		token = line_elements[i+1]
		if token[0] == "*":
			break
		if is_valid_value_param_string(token):
			(label, value) = token.split('=')
		else:
			label, value = None, token
		assigned = False
		for t in params:
			if (label is None and t['pos'] == i) or label == t['label']:
				an.update({t['dest']:convert(value, t['type'])})
				assigned = True
				break
		if assigned:
			params.pop(params.index(t))
			continue
		else:
			raise NetlistParseError, "Unknown .%s parameter: pos %d (%s=)%s" % \
			                         (atype.upper(), i, label, value)

	missing = []
	for t in params:
		if t['needed']:
			missing.append(t['label'])
	if len(missing):
		raise NetlistParseError, \
		      "Required parameters are missing: %s" % ("".join(line_elements))
	# load defaults for unsupplied parameters
	for t in params:
		an.update({t['dest']:t['default']})

	# ad-hoc code for tran ... :(
	if an['type'] == 'tran':
		uic = int(an.pop('uic'))
		if uic == 0:
			an['x0'] = None
		elif uic == 1:
			an['x0'] = 'op'
		elif uic == 2:
			an['x0'] = 'op+ic'
		elif uic == 3:
			pass # already set by ic_label
		else:
			raise NetlistParseError, "Unknown UIC value: %d" % uic
			
	return an

def is_valid_value_param_string(astr):
	"""Has the string a form like <param_name>=<value>?
	No spaces.
	Returns: a boolean
	"""
	work_astr = astr.strip()
	if work_astr.count("=") == 1:
		ret_value = True
	else:
		ret_value = False
	return ret_value
		
def convert(astr, rtype, raise_exception=False):
	if rtype == float:
		try:
			ret = convert_units(astr)
		except ValueError, msg:
			if raise_exception:
				raise ValueError, msg
			else:
				ret = astr
	elif rtype == str:
		ret = astr
	elif rtype == bool:
		ret = convert_boolean(astr)
	elif raise_exception:
		raise ValueError, "Unknown type %s" % rtype
	else:
		ret = astr
	return ret

def parse_param_value_from_string(astr, rtype=float, raise_exception=False):
	""" Searches the string for a <param>=<value> couple and returns a list.
	if rtype is float (type), default value, the method will attempt converting 
	the value to float. If the conversion fails, a string is returned. 
	If set to str (type), a string will be returned, as if the conversion failed.

	This prevents value being '0' (str) -> converted to 0.0 -> converted (somewhere 
	else) to '0.0' (str), ending up being a new node instead fo the reference.
	
	Notice that in <param>=<value> there is no space before or after the equal sign.
	
	Returns: [param, value] where param and value are both strings.
	"""
	if not is_valid_value_param_string(astr):
		return (astr, "")
	p, v = astr.strip().split("=")
	v = convert(v, rtype, raise_exception=False)
	return p, v 
	
class NetlistParseError(Exception):
	"""Netlist parsing exception."""
	pass

def convert_boolean(value):
	"""Converts the following strings to a boolean:
	yes, 1, true to True
	no, false, 0 to False
	
	raises NetlistParserException
	
	Returns: boolean
	"""
	if value == 'no' or value == 'false' or value == '0' or value == 0:
		return_value = False
	elif value == 'yes' or value == 'true' or value == '1' or value == 1:
		return_value = True
	else:
		raise NetlistParseError, "invalid boolean: "+value
	
	return return_value

def parse_ic_directive(line, line_elements=None):
	"""Parses an ic directive and assembles a dictionary accordingly.
	"""
	if line_elements is None:
		line_elements = line.split()
	
	ic_dict = {}
	name = None
	for token in line_elements[1:]:
		if token[0] == "*":
			break
		
		(label, value) = parse_param_value_from_string(token)
		if label == "name" and name is None:
			name = value
			continue
		# the user should have specified either something like:
		# V(node)=10u
		# or something like:
		# I(Vtest)=100e-6
		ic_dict.update({label:convert_units(value)})	
		# We may decide to check if the node exists and/or if the syntax
		# is correct and raise NetlistParseError if needed.
	
	if name is None:
		raise NetlistParseError, "The 'name' parameter is missing"

	return {name:ic_dict}

def parse_sub_declaration(subckt_lines):
	"""Returns a circuit.subckt instance that holds the subckt 
	information, ready to be instantiated/called.
	"""
	index = 0
	netlist_lines = []
	connected_nodes_list = []
	for line, line_n in subckt_lines:
		if index == 0:
			line_elements = line.split()
			if line_elements[0] != '.subckt':
				raise RuntimeError, "BUG? parse_sub_declaration() \
				called on non-subckt text. (line"+str(line_n)+")"
			name = line_elements[1]
			for node_name in line_elements[2:]:
				if node_name[0] == '0':
					raise NetlistParseError, "subckt " + name + \
					" has a connection node named '0' (line"+str(line_n)+")"
				if node_name[0] == '*':
					break
				else:
					connected_nodes_list = connected_nodes_list + [node_name]
		else:
			netlist_lines = netlist_lines + [(line, "")]
		index = index + 1
	subck_inst = circuit.subckt(name, netlist_lines, connected_nodes_list)
	return subck_inst

def parse_sub_instance(line, circ, subckts_dict, line_elements=None, models=None):
	"""Parses a subckt call/instance.
	
	1. Gets name and nodes connections
	2. Looks in subckts_dict for a matching subckts_dict[name]
	3. Builds a circuit wrapper
	4. Calls main_netlist_parser() on the subcircuit code 
	   (with the wrapped circuit)
	
	Returns: a elements list
	
	"""
	if line_elements is None:
		line_elements = line.split()

	if (len(line_elements) < 2):
		raise NetlistParseError, ""
	
	param_value_dict = {}
	name = None
	
	for index in range(1, len(line_elements)):
		if line_elements[index][0] == '*':
			break
		
		(param, value) = parse_param_value_from_string(line_elements[index], rtype=str)
		param_value_dict.update({param:value})
	
	if not param_value_dict.has_key("name"):
		raise NetlistParseError, "missing 'name' in subckt call"
	if not subckts_dict.has_key(param_value_dict['name']):
		raise NetlistParseError, "subckt " + param_value_dict['name'] + " is unknown"
	
	name = param_value_dict['name']
	subckt = subckts_dict[name]
	connection_nodes_dict = {}
	
	for param, value in param_value_dict.iteritems():
		if param == 'name':
			continue
		if param in subckt.connected_nodes_list:
			connection_nodes_dict.update({param:value})
		else:
			raise NetlistParseError, "unknown node " + param
	
	# check all nodes are connected
	for node in subckt.connected_nodes_list:
		if not connection_nodes_dict.has_key(node):
			raise NetlistParseError, "unconnected subckt node " + node
	
	wrapped_circ = circuit.circuit_wrapper(circ, connection_nodes_dict, subckt.name, line_elements[0])
	
	elements_list = main_netlist_parser(wrapped_circ, subckt.code, subckts_dict, models)
	
	# Every subckt adds elemets with the _same description_ (elem.descr)
	# We modify it so that each descr is uniq for every instance
	for element in elements_list:
		element.descr = "-" + wrapped_circ.prefix + element.descr
	
	return elements_list

def parse_include_directive(line, line_elements=None):
	""".include <filename> [*comments]
	"""
	if line_elements is None:
		line_elements = line.split()

	if not len(line_elements) > 1 or \
	(len(line_elements) > 2 and not line_elements[2][0] == '*'):
		raise NetlistParseError, ""
	
	path = line_elements[1]
	if not utilities.check_file(path):
		raise RuntimeError, ""
	
	fnew = open(path, "r")
	
	return [None, path, True]
	
def join_lines(fp, line):
	"""Read the lines coming up in the file. Each line that starts with '+' is added to the 
	previous line (line continuation rule). When a line not starting with '+' is found, the 
	file is rolled back and the line is returned.
	"""
	while True:
		last_pos = fp.tell()
		next = fp.readline()
		next = next.strip().lower()
		if not next:
			break
		elif next[0] == '+':
			line += ' ' + next[1:]
		else:
			fp.seek(last_pos)
			break
	return line

