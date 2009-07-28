# -*- coding: iso-8859-1 -*-
# mosq.py
# Square law mos model
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
	mos_type: may be 'n' or 'p', identifies the mos type"""

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
		if self.mos_type == 'n':
			if  ports_v[0] - ports_v[1] >= 0:
				v1 = ports_v[0]
				v2 = ports_v[1]
				reversed_factor = 1
			else:
				v1 = ports_v[1]
				v2 = ports_v[0]
				reversed_factor = -1
			sign = +1
	
		if self.mos_type == 'p':
			sign = -1
			if  ports_v[0] - ports_v[1] <= 0:
				v1 = sign * ports_v[0]
				v2 = sign * ports_v[1]
				reversed_factor = 1
			else:
				v1 = sign * ports_v[1]
				v2 = sign * ports_v[0]
				reversed_factor = -1
			
		return reversed_factor*sign*self._get_id(vgs=v1, vgd=v2)

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
		if self.mos_type == 'n':
			if ports_v[0] - ports_v[1] >= 0:
				v1 = ports_v[0]
				v2 = ports_v[1]
				reversed_factor = +1
			else:
				v1 = ports_v[1]
				v2 = ports_v[0]
				if port_index == 1:
					port_index = 0
				elif port_index == 0:
					port_index = 1
				reversed_factor = -1	
		elif self.mos_type == 'p':
			sign = -1
			if  ports_v[0] - ports_v[1] <= 0:
				v1 = sign * ports_v[0]
				v2 = sign * ports_v[1]
				reversed_factor = +1
			else:
				v1 = sign * ports_v[1]
				v2 = sign * ports_v[0]
				reversed_factor = -1
				if port_index == 1:
					port_index = 0
				elif port_index == 0:
					port_index = 1
		
		return reversed_factor*self._get_g(port_index, vgs=v1, vgd=v2)

	def _get_g(self, port_index, vgs, vgd):
		if port_index == 0: #vgs
			if vgs <= self.vt and vgd <= self.vt:
				# no channel at both sides
				gdr = 0
			elif vgs > self.vt and vgd > self.vt:
				# zona ohmica: channel at both sides
				gdr =  self.kp * (self.w / self.l) * (vgs - self.vt)
			elif vgs > self.vt and vgd <= self.vt:
				# zona di saturazione: canale al s
				gdr =  self.kp * (vgs - self.vt) * (self.w / self.l) * (1 - self.lambd*(vgd - self.vt))
			else:
				# zona di saturazione: canale al d # era: 3*vgs - vds - 3*self.vt
				gdr =  0.5 * self.kp * ((vgd - self.vt)**2) * (self.w / self.l) * self.lambd*vgs
		elif port_index == 1: #vgd
			if vgs <= self.vt and vgd <= self.vt:
				# no channel at both sides
				gdr = 0
			elif vgs > self.vt and vgd > self.vt:
				# zona ohmica: channel at both sides
				gdr = self.kp * (self.vt - vgd) * (self.w / self.l)
			elif vgs > self.vt and vgd <= self.vt:
				# zona di saturazione: canale al s
				gdr = - 0.5 * self.kp * ((vgs - self.vt)**2) * (self.w / self.l) * self.lambd
			else:
				# zona di saturazione: canale al d
				gdr =  - self.kp * (vgd - self.vt) * (self.w / self.l) * (1 - self.lambd*(vgs - self.vt))
		
		return gdr

	def _get_id(self, vgs, vgd):
		if vgs <= self.vt and vgd <= self.vt:
			# no channel at both sides
			idr = 0
			if self._debug:
				print "M"+self.descr+":", "vgs:", str(vgs), "vgd:", str(vgd), "vds:", str(vgs-vgd), "OFF"
		elif vgs > self.vt and vgd > self.vt:
			# zona ohmica: channel at both sides
			idr = .5 * self.kp * (vgs - vgd)*(-2*self.vt + vgs + vgd) * (self.w / self.l)
			if self._debug:
				print "M"+self.descr+":", "vgs:", str(vgs), "vgd:", str(vgd), "vds:", str(vgs-vgd), "OHM", "idr:", str(idr)
		elif vgs > self.vt and vgd <= self.vt:
			# zona di saturazione: canale al s
			idr =  0.5 * self.kp * ((vgs - self.vt)**2) * (self.w / self.l) * (1 - self.lambd*(vgd - self.vt))
			if self._debug:
				print "M"+self.descr+":", "vgs:", str(vgs), "vgd:", str(vgd), "vds:", str(vgs-vgd), "SAT", "idr:", str(idr)
		else:
			# zona di saturazione: canale al d
			idr =  -0.5 * self.kp * ((vgd - self.vt)**2) * (self.w / self.l) * (1 - self.lambd*(vgs - self.vt))
			if self._debug:
				print "M"+self.descr+":", "vgs:", str(vgs), "vgd:", str(vgd), "vds:", str(vgs-vgd), "SAT", "idr:", str(idr)
		return idr
	
#	def print_op_info(self, )

if __name__ == '__main__': 
	import numpy
	mymos = mosq(nd=2, ng=1, ns=0, kp=1, w=1, l=1, vt=1, lambd=0.5, mos_type='p')
	mymos.descr = "TEST"
	vgs_v = (numpy.array(range(400))-200)/100.0
	vgd_v = [-2, -1.5, -1, -0.5, -0.2, 0, 0.2, 0.5, 1, 1.5, 2]
	for vgd in vgd_v:
		fp = open("mos_test-vgs_sweep-vgd_"+str(vgd)+".txt", 'w')
		fp.write("#VGS\tVGD\tI\tgm\tgos\n")
		for vgs in vgs_v:
			fp.write(str(vgs)+"\t"+str(vgd)+"\t"+str(mymos.i([vgs, vgd]))+"\t"+str(mymos.g([vgs, vgd], 0))+"\t"+str(mymos.g([vgs, vgd], 1))+"\n")
		fp.close()

