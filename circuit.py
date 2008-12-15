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
Holds the circuit class and various element classes.

General form of a _non linear_ element class:

The class must provide:

1. Element terminals:
elem.n1 # the anode of the element
elem.n2 # the cathode of the element

Notice: a positive current is a current that flows into the anode and out of
the cathode. This convention is used throughout the simulator.

2. elem.get_ports()
This method must return a tuple of pairs of nodes. Something like:
((na, nb), (nc, nd), (ne, nf), ... )

Each pair of nodes is used to determine a voltage that has effect on the
current. I referred to them as a 'port' though it may not be a good idea.

For example, an nmos has:
((n_gate, n_source), (n_drain, n_source))

The positive terminal is the first.

From that, the calling method builds a voltage vector corresponding to the
ports vector:
voltages_vector = ( Va-Vb, Vc-Vd, Ve-Vf, ...)

That's passed to:
3. elem.i(voltages_vector, time)

It returns the current flowing into the element if the voltages specified in
the voltages_vector are applied to its ports, at the time given.

4. elem.g(voltages_vector, port_index, time) is similar, but returns the 
differential transconductance between the port at position port_index in the 
ports_vector (see 2) and the current, when the operating point is specified by
the voltages in the voltages_vector. 

5. elem.is_nonlinear
A non linear element must have a elem.is_nonlinear field set to True.

Recommended:
1. A non linear element may have a list/tuple of the same length of its 
ports_vector in which there are the recommended guesses for dc analysis. 
Eg Vgs is set to Vt in mosfets.
This is obviously useless for linear devices.

2. Every element should have a meaningful __str__ method. 
It must return a line of paramaters without n1 n2, because a element cannot 
know the external names of its nodes. It is used to print the parsed netlist.

"""

import math
import constants, options

# will be added here by netlist_parser
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
	See netlist_parser for examples of how I did it.

	The code to add a element should be something like:

		anode = my_circuit.add_node_to_circuit(ext_name_of_anode)
		cathode = my_circuit.add_node_to_circuit(ext_name_of_cathode)
		resistance = 1e3 # 1Kohm
		# new_element = circuit.element( ... ) 
		new_element = circuit.resitor(n1=anode, n2=cathode, R=resistance)
		my_circuit.elements.append(new_element)
	
	3. Internal only nodes
	The number of internal only nodes (added automatically by the simulator)
	is hold in my_circuit.internal_nodes
	That value shouldn't be changed by hand.
	"""
	title = ""
	nodes_dict = {} # {int_node:ext_node}
	#_reverse_dict = {}
	elements = []
	internal_nodes = 0
	
	
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
			ports = elem.get_ports()
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
	
## NOTES ON ELEM CLASSES
# see above ^^^^^

#class generic:
	#"""This class is for debugging purposes only."""
	##generic element
	#is_nonlinear = False
	#n1 = None
	#n2 = None
	#def __init__(self, n1, n2, is_nonlinear):
		#self.n1 = n1
		#self.n2 = n2
		#self.is_nonlinear = is_nonlinear
	##must be called to define the element!
	#def set_char(self, i_function=None, g_function=None):
		#if i_function:
			#self.i = i_function
		#if g_function:
			#self.g = g_function
	##def g(self, v):
	##	return 1/self.R
	##def i(self, v):
	##	return 0

class resistor:
	letter_id = "r"
	is_nonlinear = False
	def __init__(self, n1, n2, R):
		self.R = R
		self.n1 = n1
		self.n2 = n2
	def __str__(self):
		return str(self.R)
	def g(self, v, time=0):
		return 1.0/self.R
	def i(self, v, time=0):
		return 0
class capacitor:
	letter_id = "c"
	is_nonlinear = False
	def __init__(self, n1, n2, C, ic=None):
		self.C = C
		self.n1 = n1
		self.n2 = n2
		self.ic = ic
	def __str__(self):
		return str(self.C)
	def g(self, v, time=0):
		return 0
	def i(self, v, time=0):
		return 0
	def d(self, v, time=0):
		return self.C
class inductor:
	letter_id = "l"
	is_nonlinear = False
	def __init__(self, n1, n2, L, ic=None):
		self.L = L
		self.n1 = n1
		self.n2 = n2
		self.ic = ic
	def __str__(self):
		return str(self.L)


#########################
## NON LINEAR ELEMENTS
#########################

class diode:
	letter_id = "d"
	is_nonlinear = True
	dc_guess = [0.4]
	def __init__(self, n1, n2, Io=1.0e-14, m=1, T=None, ic=None):
		if Io:
			self.Io = Io
		else:
			self.Io = 1.0e-14
		
		self.n1 = n1
		self.n2 = n2
		self.ports = ((self.n1, self.n2),)
		if m:
			self.m = m
		else:
			self.m = 1
		if T:
			self.T = T
		else:
			self.T = constants.T
		
		if ic is not None: #fixme
			print "(W): ic is ignored in diodes."
		self.ic = ic #per ora inutilizzato
	def __str__(self):
		rep = "Is="+str(self.Io)+" m="+str(self.m)+" T="+str(self.T)
		if self.ic is not None:
			rep = rep + " ic="+str(self.ic)
		return rep
	def get_ports(self):
		return self.ports
	def i(self, ports_v, time=0): #with gmin added
		v = ports_v[0]
		return self.Io*(math.exp((constants.e*v)/(self.m*constants.k*self.T))-1) + options.gmin*v
	def g(self, ports_v, port_index, time=0):
		if not port_index == 0: 
			raise Exception, "Attepted to evaluate a diode's gm on a unknown port."
		return (self.i(ports_v)*constants.e)/(self.m*constants.k*self.T) + options.gmin

class mosq:
	"""Square law MOS model
	
	nd: drain node
	ng: gate node
	ns: source node
	kp: kp = u0 * C'ox
	w: channel width in meters
	l: channel length in meters
	vt: threshold voltage
	lambd: lambd = 1/Va channel length modulation
	mos_type: may be 'n' or 'p', identifies the mos type
	"""
	letter_id = "m"
	is_nonlinear = True
	dc_guess = None #defined in init
	descr = None
	
	_debug = False
	
	def __init__(self, nd, ng, ns, kp, w, l, vt, lambd=0, mos_type='n'):
		self.ng = ng
		self.n1 = nd
		self.n2 = ns
		self.kp = kp
		self.w = float(w)
		self.l = float(l)
		self.vt = vt
		self.ports = ((self.ng, self.n2), (self.ng, self.n1))
		if mos_type != 'n' and mos_type != 'p':
			raise Exception, "Unknown mos type: " + str(mos_type)
		self.mos_type = mos_type.lower()
		self.lambd = lambd * w / l
		self.dc_guess = [self.vt*(1.1)*(1-2*(self.mos_type=='p')), 0]
	
	def get_ports(self):
		"""Returns a tuple of tuples of ports nodes, as:
		(port0, port1, port2...)
		Where each port is in the form:
		port0 = (nplus, nminus)
		"""
		return self.ports
	
	def __str__(self):
		rep = "type=" + self.mos_type + " kp=" + str(self.kp)+ " vt="+ str(self.vt)+ " w="+ str(self.w)+ " l=" + str(self.l) + \
		" lambda="+ str(self.lambd)
		return rep
	

	def i(self, ports_v, time=0):
		"""Returns the current flowing in the element with the voltages applied
		as specified in the ports_v vector.
		
		ports_v: a list in the form: [voltage_across_port0, voltage_across_port1, ...]
		time: the simulation time at which the evaluation is performed. Set it to
		None during DC analysis.
		
		"""
		v1 = ports_v[0]
		v2 = ports_v[1]
		
		if self.mos_type == 'p':
			v1 = -v1
			v2 = -v2
			sign = -1
		else:
			sign = +1
			
			
		if v1 <= self.vt and v2 <= self.vt:
			# no channel at both sides
			idr = 0
			if self._debug:
				print "M"+self.descr+":", "vgs:", str(v1), "vgd:", str(v2), "vds:", str(v1-v2), "OFF"
		elif v1 > self.vt and v2 > self.vt:
			# zona ohmica: channel at both sides
			idr = .5 * self.kp * (v1 - v2)*(-2*self.vt + v1 + v2) * (self.w / self.l)
			if self._debug:
				print "M"+self.descr+":", "vgs:", str(v1), "vgd:", str(v2), "vds:", str(v1-v2), "OHM", "idr:", str(idr)
		elif v1 > self.vt and v2 <= self.vt:
			# zona di saturazione: canale al s
			idr =  0.5 * self.kp * ((v1 - self.vt)**2) * (self.w / self.l) * (1 - self.lambd*(v2 - self.vt))
			if self._debug:
				print "M"+self.descr+":", "vgs:", str(v1), "vgd:", str(v2), "vds:", str(v1-v2), "SAT", "idr:", str(idr)
		else:
			# zona di saturazione: canale al d
			idr =  -0.5 * self.kp * ((v2 - self.vt)**2) * (self.w / self.l) * (1 - self.lambd*(v1 - self.vt))
			if self._debug:
				print "M"+self.descr+":", "vgs:", str(v1), "vgd:", str(v2), "vds:", str(v1-v2), "SAT", "idr:", str(idr)

		return sign*idr
	
	def g(self, ports_v, port_index, time=0):
		"""
		Returns the differential (trans)conductance rs the port specified by port_index
		when the element has the voltages specified in ports_v across its ports,
		at (simulation) time.
		
		ports_v: a list in the form: [voltage_across_port0, voltage_across_port1, ...]
		port_index: an integer, 0 <= port_index < len(self.get_ports())
		time: the simulation time at which the evaluation is performed. Set it to
		None during DC analysis.
		"""
		v1 = ports_v[0]
		v2 = ports_v[1]
		vds = v1 - v2
		if self.mos_type == 'p':
			v1 = -v1
			v2 = -v2
			
		if port_index == 0: #vgs
			if v1 <= self.vt and v2 <= self.vt:
				# no channel at both sides
				gdr = 0
			elif v1 > self.vt and v2 > self.vt:
				# zona ohmica: channel at both sides
				gdr =  self.kp * (self.w / self.l) * (v1 - self.vt)
			elif v1 > self.vt and v2 <= self.vt:
				# zona di saturazione: canale al s
				gdr =  self.kp * (v1 - self.vt) * (self.w / self.l) * (1 - self.lambd*(v2 - self.vt))
			else:
				# zona di saturazione: canale al d # era: 3*vgs - vds - 3*self.vt
				gdr =  0.5 * self.kp * ((v2 - self.vt)**2) * (self.w / self.l) * self.lambd*v1
		elif port_index == 1: #vgd
			if v1 <= self.vt and v2 <= self.vt:
				# no channel at both sides
				gdr = 0
			elif v1 > self.vt and v2 > self.vt:
				# zona ohmica: channel at both sides
				gdr = self.kp * (self.vt - v2) * (self.w / self.l)
			elif v1 > self.vt and v2 <= self.vt:
				# zona di saturazione: canale al s
				gdr = - 0.5 * self.kp * ((v1 - self.vt)**2) * (self.w / self.l) * self.lambd
			else:
				# zona di saturazione: canale al d
				gdr =  - self.kp * (v2 - self.vt) * (self.w / self.l) * (1 - self.lambd*(v1 - self.vt))
		
		return gdr



################
##  SOURCES
################

class isource:
	"""Generic (ideal) current source:
	Defaults to a DC current source. To implement a time-varying source:
	set _time_function to an appropriate function(time) and is_timedependent=True
	
	n1: + node
	n2: - node
	idc: DC current (A)
	
	Note: if DC voltage is set and is_timedependent == True, idc will be returned
	if the current is evaluated in a DC analysis. This may be useful to simulate a OP
	and then perform a transient analysis with the OP as starting point.
	Otherwise the value in t=0 is used for DC analysis.
	"""
	letter_id = "i"
	is_nonlinear = False
	is_timedependent = False
	_time_function = None
	def __init__(self, n1, n2, idc=None):
		self.idc = idc
		self.n1 = n1
		self.n2 = n2
		
	def __str__(self):
		rep = ""
		if self.idc is not None:
			rep = rep + "type=idc idc="+str(self.idc) + " "
		if self.is_timedependent:
			rep = rep + str(self._time_function)
		return rep

	def I(self, time=None):
		"""Returns the current in A at the time supplied.
		If time is not supplied, or set to None, or the source is DC, returns idc
		
		This simulator uses Normal convention: 
		A positive currents flows in a element from the + node to the - node
		"""
		if not self.is_timedependent or (self._time_function == None) or (time==None and self.idc is not None):
			return self.idc
		else:
			return self._time_function.value(time)
class vsource:
	"""Generic (ideal) voltage source:
	Defaults to a DC voltage source. To implement a time-varying source:
	set _time_function to an appropriate function(time) and is_timedependent=True
	
	n1: + node
	n2: - node
	vdc: DC voltage (V)
	
	Note: if DC voltage is set and is_timedependent == True, vdc will be returned
	if the voltage is evaluated in a DC analysis. This may be useful to simulate a OP
	and then perform a transient analysis with the OP as starting point.
	Otherwise the value in t=0 is used for DC analysis.
	"""
	letter_id = "v"
	is_nonlinear = False
	is_timedependent = False
	_time_function = None
	dc_guess = None #defined in init
	def __init__(self, n1, n2, vdc=None):
		self.vdc = vdc
		self.n1 = n1
		self.n2 = n2
		if vdc is not None:
			self.dc_guess = [self.vdc]
	
	def __str__(self):
		rep = ""
		if self.vdc is not None:
			rep = rep + "type=vdc vdc="+str(self.vdc) + " "
		if self.is_timedependent:
			rep = rep + str(self._time_function)
		return rep
	
	def V(self, time=None):
		"""Returns the voltage in V at the time supplied.
		If time is not supplied, or set to None, or the source is DC, returns vdc"""
		if not self.is_timedependent or \
		(self._time_function is None) or \
		(time is None and self.vdc is not None):
			return self.vdc
		else:
			return self._time_function.value(time)

class evsource:
	"""Linear voltage controlled voltage source (ideal)

	Source port is a open circuit, dest. port is a ideal voltage source:
	(Vn1 - Vn2) = alpha * (Vsn1 - Vsn2)
	
	n1: + node, output port
	n2: - node, output port
	sn1: + node, source port
	sn2: - node, source port
	alpha: prop constant between voltages
	
	"""
	letter_id = "e"
	is_nonlinear = False
	def __init__(self, n1, n2, sn1, sn2, alpha):
		self.alpha = alpha
		self.n1 = n1
		self.n2 = n2
		self.sn1 = sn1
		self.sn2 = sn2
	def __str__(self):
		return "alpha="+str(self.alpha)


class gisource:
	letter_id = "g"
	"""Linear voltage controlled current source
	
	Source port is a short circuit, dest. port is a ideal current source:
	Io = alpha * Is
		
	Where a positive I enters in n+ and exits from n-
	
	n1: + node, output port
	n2: - node, output port
	sn1: + node, source port
	sn2: - node, source port
	alpha: prop constant between currents
	
	"""
	is_nonlinear = False
	def __init__(self, n1, n2, sn1, sn2, alpha):
		self.alpha = alpha
		self.n1 = n1
		self.n2 = n2
		self.sn1 = sn1
		self.sn2 = sn2
	def __str__(self):
		return "alpha="+str(self.alpha)
	
class hvsource: #fixme
	"""Linear current controlled voltage source
	fixme: todo
	
	"""
	letter_id = "h"
	is_nonlinear = False
	def __init__(self, n1, n2, sn1, sn2, alpha):
		print "hvsource not implemented. BUG"
		self.alpha = alpha
		self.n1 = n1
		self.n2 = n2
		self.sn1 = sn1
		self.sn2 = sn2
		import sys
		sys.exit(1)
	def __str__(self):
		raise Exception, "hvsource not implemented. BUG"

# NEEDS TO BE CALLED hvsource, or search for it and modify appropriately



#########################################
# Functions for time dependent sources  #
#########################################

class pulse:
	#PULSE(V1 V2 TD TR TF PW PER)
	_type = "V"
	v1 = None
	v2 = None
	td = None
	per = None
	tr = None
	tf = None
	pw = None
	def __init__(self, v1=None, v2=None, td=None, tr=None, pw=None, tf=None, per=None):
		self.v1 = v1
		self.v2 = v2
		self.td = td
		self.per = per
		self.tr = tr
		self.tf = tf
		self.pw = pw
	def value(self, time):
		if not self.ready():
			print "Error: pulse function not well defined. This is a bug."
		time = time - self.per*int(time/self.per)
		if time < self.td:
			return self.v1
		elif time < self.td+self.tr:
			return self.v1 + ((self.v2-self.v1)/(self.tr))*(time - self.td)
		elif time < self.td+self.tr+self.pw:
			return self.v2
		elif time < self.td+self.tr+self.pw+self.tf:
			return self.v2 + ((self.v1-self.v2)/(self.tf))*(time - (self.td+self.tr+self.pw))
		else:
			return self.v1
	def ready(self):
		if self.v1 == None or self.v2 == None or self.td == None or self.tr == None or self.pw == None or \
		self.tf == None or self.per == None:
			return False
		else:
			return True
	def __str__(self):
		return "type=pulse " + \
		self._type.lower() + "1="+str(self.v1) + " " + \
		self._type.lower() + "2=" + str(self.v2) + \
		" td=" + str(self.td) + " per=" + str(self.per) + \
		" tr=" + str(self.tr) + " tf=" + str(self.tf) + \
		" pw=" + str(self.pw)
class sin:
	#SIN(VO VA FREQ TD THETA)
	td = None
	vo = None
	va  = None
	freq  = None
	theta = None
	_type = "V"
	def __init__(self, vo=None, va=None, freq=None, td=None, theta=None):
		self.vo = vo
		self.va = va
		self.freq = freq
		self.td = td
		self.theta = theta
	def value(self, time):
		if not self.ready():
			print "Error: sin function not well defined. This is a bug."
		if time < self.td:
			return self.vo
		elif self.theta:
			return self.vo + self.va * math.exp(-1*(time-self.td)/self.theta) * math.sin(2*math.pi*self.freq*(time+self.td))
		else:
			return self.vo + self.va * math.sin(2*math.pi*self.freq*(time+self.td))
	def ready(self):
		if self.vo == None or self.va == None or self.freq == None or self.td == None or self.theta == None:
			return False
		else:
			return True
	def __str__(self):
		return "type=sin " + \
		self._type.lower() + "o=" + str(self.vo) + " " + \
		self._type.lower() +"a=" + str(self.va) + \
		" freq=" + str(self.freq) + " theta=" + str(self.theta) + \
		" td=" + str(self.td)
	

class exp:
	# EXP(V1 V2 TD1 TAU1 TD2 TAU2)
	v1 = None
	v2 = None
	td1  = None
	tau1  = None
	td2 = None
	tau2 = None
	_type = "V"
	def __init__(self, v1=None, v2=None, td1=None, tau1=None, td2=None, tau2=None):
		self.v1 = v1
		self.v2 = v2
		self.td1 = td1
		self.tau1 = tau1
		self.td2 = td2
		self.tau2 = tau2
	def value(self, time):
		if not self.ready():
			print "Error: exp function not well defined. This is a bug."
		if time < self.td1:
			return self.v1
		elif time < self.td2:
			return self.v1+(self.v2 - self.v1) * (1-math.exp(-1 * (time - self.td1)/self.tau1))
		else:
			return self.v1 + (self.v2 - self.v1) * (1 - math.exp(-1 * (time - self.td1 ) / self.tau1))+(self.v1 - self.v2 ) * ( 1 - math.exp(-1 * (time - self.td2) / self.tau2))
	def ready(self):
		if self.v1 == None or self.v2 == None or self.td1 == None or self.tau1 == None or self.td2 == None \
		or self.tau2 == None:
			return False
		return True
	def __str__(self):
		return "type=exp " + \
		self._type.lower() + "1=" + str(self.v1) + " " + \
		self._type.lower() + "2=" + str(self.v2) + \
		" td1="+str(self.td1) + " td2=" + str(self.td2) + \
		" tau1=" + str(self.tau1) + " tau2=" + str(self.tau2)
														

# STATIC METHODS
def is_elem_voltage_defined(elem):
	"""Returns: 
	True se the elem is a vsource, inductor, evsource or hvsource
	False otherwise.
	"""
	if isinstance(elem, vsource) or isinstance(elem, evsource) or \
	isinstance(elem, hvsource) or isinstance(elem, inductor):
		return True
	else:
		return False

class NodeNotFoundError(Exception):
	"""Circuit Node exception."""
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
	circ = None
	subckt_node_filter_dict = {}
	prefix = ""
	def __init__(self, circ, connection_nodes_dict, subckt_name, subckt_label):
		self.circ = circ
		self.prefix = subckt_label + "-" + subckt_name + "-"
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
