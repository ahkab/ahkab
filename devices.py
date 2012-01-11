# -*- coding: iso-8859-1 -*-
# devices.py
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

import math
import constants, printing

"""
Contains various element classes.

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
	def get_op_info(self, ports_v):
		vn1n2 = float(ports_v[0][0])
		in1n2 = float(ports_v[0][0]/self.R)
		power = float(ports_v[0][0]**2/self.R)
		arr = [[self.letter_id.upper()+self.descr,"V(n1-n2):", vn1n2, "[V]", "I(n2-n1):", in1n2, "[A]", "P:", power, "[W]"]]
		strarr = printing.table_setup(arr)
		return strarr
	def print_op_info(self, ports_v):
		print self.get_op_info(ports_v),
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
	def get_op_info(self, ports_v):
		vn1n2 = float(ports_v[0][0])
		qn1n2 = float(ports_v[0][0]*self.C)
		energy = float(.5*ports_v[0][0]**2*self.C)
		arr = [[self.letter_id.upper()+self.descr,"V(n1-n2):", vn1n2, "[V]", "Q:", qn1n2, "[C]", "E:", energy, "[J]"]]
		strarr = printing.table_setup(arr)
		return strarr
	def print_op_info(self, ports_v):
		print self.get_op_info(ports_v),
class inductor:
	letter_id = "l"
	is_nonlinear = False
	def __init__(self, n1, n2, L, ic=None):
		self.L = L
		self.n1 = n1
		self.n2 = n2
		self.ic = ic
		self.coupling_devices = []
	def __str__(self):
		return str(self.L)

class inductor_coupling:
	letter_id = "k"
	is_nonlinear = False
	def __init__(self, L1, L2, K, M):
		self.K = K
		self.M = M
		self.L1 = L1
		self.L2 = L2
	def __str__(self):
		return "%s %s %g" % (self.L1, self.L2, self.K)
	def get_other_inductor(self, Lselected):
		Lret = None
		if Lselected.upper() == self.L1.upper():
			Lret = self.L2
		elif Lselected.upper() == self.L2.upper():
			Lret = self.L1
		if Lret is None:
			raise Exception, "Mutual inductors bug."
		return Lret

#########################
## NON LINEAR ELEMENTS
#########################

class diode:
	letter_id = "d"
	is_nonlinear = True
	dc_guess = [0.4]
	def __init__(self, n1, n2, Io=1.0e-14, m=1, T=None, ic=None, off=False):
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
			self.T = None
		
		if ic is not None: #fixme
			print "(W): ic support in diodes is very experimental."
			self.dc_guess = ic
		self.ic = ic
		self.off = off
		if self.off:
			if self.ic is None:
				self.ic = 0
			else:
				print "(W): IC statement in diodes takes precedence over OFF."
				print "(W): If you are performing a transient simulation with uic=2,"
				print "(W): you may want to check the initial value."

	def _get_T(self):
		if self.T is None: 
			return constants.T
		else:
			return self.T

	def __str__(self):
		T = self._get_T()
		rep = "Is="+str(self.Io)+" m="+str(self.m)+" T="+str(T)
		if self.ic is not None:
			rep = rep + " ic="+str(self.ic)
		return rep

	def get_output_ports(self):
		return self.ports

	def get_drive_ports(self, op):
		return self.ports

	def i(self, op_index, ports_v, time=0): #with gmin added
		v = ports_v[0]
		return self.Io*(math.exp(v/(self.m*constants.Vth(self._get_T())))-1)

	def g(self, op_index, ports_v, port_index, time=0):
		if not port_index == 0: 
			raise Exception, "Attepted to evaluate a diode's gm on a unknown port."
		return (self.i(op_index, ports_v)/constants.Vth(self._get_T())/self.m)

	def get_op_info(self, ports_v_v):
		vn1n2 = float(ports_v_v[0][0])
		idiode = self.i(0, (vn1n2,))
		gmdiode = self.g(0, (vn1n2,), 0)
		info = ["V(n1-n2): %g ", vn1n2, "[V]", "I(n1-n2):", idiode, "[A]", "P:", vn1n2*idiode, "g:", gmdiode, "[A/V]", "T:", self._get_T(), "K" ]
		arr = [[self.letter_id.upper()+self.descr] + info]
		strarr = printing.table_setup(arr)
		return strarr
	def print_op_info(self, ports_v):
		print self.get_op_info(ports_v),

#class mosq:
#	"""Square law MOS model
#	
#	nd: drain node
#	ng: gate node
#	ns: source node
#	kp: kp = u0 * C'ox
#	w: channel width in meters
#	l: channel length in meters
#	vt: threshold voltage
#	lambd: lambd = 1/Va channel length modulation
#	mos_type: may be 'n' or 'p', identifies the mos type
#	"""
#	letter_id = "m"
#	is_nonlinear = True
#	dc_guess = None #defined in init
#	descr = None
#	
#	_debug = False
#	
#	def __init__(self, nd, ng, ns, kp, w, l, vt, lambd=0, mos_type='n'):
#		self.ng = ng
#		self.n1 = nd
#		self.n2 = ns
#		self.kp = kp
#		self.w = float(w)
#		self.l = float(l)
#		self.vt = vt
#		self.ports = ((self.ng, self.n2), (self.ng, self.n1))
#		if mos_type != 'n' and mos_type != 'p':
#			raise Exception, "Unknown mos type: " + str(mos_type)
#		self.mos_type = mos_type.lower()
#		self.lambd = lambd * w / l
#		self.dc_guess = [self.vt*(1.1)*(1-2*(self.mos_type=='p')), 0]
#	
#	def get_ports(self):
#		"""Returns a tuple of tuples of ports nodes, as:
#		(port0, port1, port2...)
#		Where each port is in the form:
#		port0 = (nplus, nminus)
#		"""
#		return self.ports
#	
#	def __str__(self):
#		rep = "type=" + self.mos_type + " kp=" + str(self.kp)+ " vt="+ str(self.vt)+ " w="+ str(self.w)+ " l=" + str(self.l) + \
#		" lambda="+ str(self.lambd)
#		return rep
#	
#
#	def i(self, ports_v, time=0):
#		"""Returns the current flowing in the element with the voltages applied
#		as specified in the ports_v vector.
#		
#		ports_v: a list in the form: [voltage_across_port0, voltage_across_port1, ...]
#		time: the simulation time at which the evaluation is performed. Set it to
#		None during DC analysis.
#		
#		"""
#		v1 = ports_v[0]
#		v2 = ports_v[1]
#		
#		if self.mos_type == 'p':
#			v1 = -v1
#			v2 = -v2
#			sign = -1
#		else:
#			sign = +1
#			
#			
#		if v1 <= self.vt and v2 <= self.vt:
#			# no channel at both sides
#			idr = 0
#			if self._debug:
#				print "M"+self.descr+":", "vgs:", str(v1), "vgd:", str(v2), "vds:", str(v1-v2), "OFF"
#		elif v1 > self.vt and v2 > self.vt:
#			# zona ohmica: channel at both sides
#			idr = .5 * self.kp * (v1 - v2)*(-2*self.vt + v1 + v2) * (self.w / self.l)
#			if self._debug:
#				print "M"+self.descr+":", "vgs:", str(v1), "vgd:", str(v2), "vds:", str(v1-v2), "OHM", "idr:", str(idr)
#		elif v1 > self.vt and v2 <= self.vt:
#			# zona di saturazione: canale al s
#			idr =  0.5 * self.kp * ((v1 - self.vt)**2) * (self.w / self.l) * (1 - self.lambd*(v2 - self.vt))
#			if self._debug:
#				print "M"+self.descr+":", "vgs:", str(v1), "vgd:", str(v2), "vds:", str(v1-v2), "SAT", "idr:", str(idr)
#		else:
#			# zona di saturazione: canale al d
#			idr =  -0.5 * self.kp * ((v2 - self.vt)**2) * (self.w / self.l) * (1 - self.lambd*(v1 - self.vt))
#			if self._debug:
#				print "M"+self.descr+":", "vgs:", str(v1), "vgd:", str(v2), "vds:", str(v1-v2), "SAT", "idr:", str(idr)
#
#		return sign*idr
#	
#	def g(self, ports_v, port_index, time=0):
#		"""
#		Returns the differential (trans)conductance rs the port specified by port_index
#		when the element has the voltages specified in ports_v across its ports,
#		at (simulation) time.
#		
#		ports_v: a list in the form: [voltage_across_port0, voltage_across_port1, ...]
#		port_index: an integer, 0 <= port_index < len(self.get_ports())
#		time: the simulation time at which the evaluation is performed. Set it to
#		None during DC analysis.
#		"""
#		v1 = ports_v[0]
#		v2 = ports_v[1]
#		vds = v1 - v2
#		if self.mos_type == 'p':
#			v1 = -v1
##			v2 = -v2
#			
#		if port_index == 0: #vgs
#			if v1 <= self.vt and v2 <= self.vt:
#				# no channel at both sides
#				gdr = 0
#			elif v1 > self.vt and v2 > self.vt:
#				# zona ohmica: channel at both sides
#				gdr =  self.kp * (self.w / self.l) * (v1 - self.vt)
#			elif v1 > self.vt and v2 <= self.vt:
#				# zona di saturazione: canale al s
#				gdr =  self.kp * (v1 - self.vt) * (self.w / self.l) * (1 - self.lambd*(v2 - self.vt))
#			else:
#				# zona di saturazione: canale al d # era: 3*vgs - vds - 3*self.vt
#				gdr =  0.5 * self.kp * ((v2 - self.vt)**2) * (self.w / self.l) * self.lambd*v1
#		elif port_index == 1: #vgd
#			if v1 <= self.vt and v2 <= self.vt:
#				# no channel at both sides
#				gdr = 0
#			elif v1 > self.vt and v2 > self.vt:
#				# zona ohmica: channel at both sides
#				gdr = self.kp * (self.vt - v2) * (self.w / self.l)
#			elif v1 > self.vt and v2 <= self.vt:
#				# zona di saturazione: canale al s
#				gdr = - 0.5 * self.kp * ((v1 - self.vt)**2) * (self.w / self.l) * self.lambd
#			else:
#				# zona di saturazione: canale al d
#				gdr =  - self.kp * (v2 - self.vt) * (self.w / self.l) * (1 - self.lambd*(v1 - self.vt))
#		
#		return gdr



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
	def __init__(self, n1, n2, idc=None, abs_ac=None, arg_ac=0):
		self.idc = idc
		self.abs_ac = abs_ac
		self.arg_ac = arg_ac
		self.n1 = n1
		self.n2 = n2
		
	def __str__(self):
		rep = ""
		if self.idc is not None:
			rep = rep + "type=idc idc="+str(self.idc) + " "
		if self.abs_ac is not None:
			rep = rep + "iac="+str(self.abs_ac) + " " + "arg="+str(self.arg_ac) + " "
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
	def __init__(self, n1, n2, vdc=None, abs_ac=None, arg_ac=0):
		self.vdc = vdc
		self.n1 = n1
		self.n2 = n2
		self.abs_ac = abs_ac
		self.arg_ac = arg_ac
		if vdc is not None:
			self.dc_guess = [self.vdc]
	
	def __str__(self):
		rep = ""
		if self.vdc is not None:
			rep = rep + "type=vdc vdc="+str(self.vdc) + " "
		if self.abs_ac is not None:
			rep = rep + "vac="+str(self.abs_ac) + " " + "arg="+str(self.arg_ac) + " "
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
			printing.print_general_error("Error: sin function not well defined. This is a bug.")
		if time < self.td:
			return self.vo
		elif self.theta:
			return self.vo + self.va * math.exp(-1*(time-self.td)/self.theta) * math.sin(2*math.pi*self.freq*(time-self.td))
		else:
			return self.vo + self.va * math.sin(2*math.pi*self.freq*(time-self.td))
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
			printing.print_general_error("Error: exp function not well defined. This is a bug.")
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
