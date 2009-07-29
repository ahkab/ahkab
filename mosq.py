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

This MOS Model follows the Square Law Mos Model:
[Vds > 0 in the following, transistor type: N]
1. No subthreshold conduction.
   Vgs < Vt
   Id = 0
2. Ohmic region of operation
   Vgs > Vt
   Vgd > Vt
   Id = k w/l ((vgs-vt)vds - vds^2/2)
3. Saturation region of operation
   Vgs > Vt
   Vgs < Vt
   Id = 1/2 k w/l (vgs-vt)^2 * (1 + lambd*(vds-vgs+vt))

Two ports: Vgs, Vgd

Positive current flows out of source and into drain.

The model is symmetric, Source and Drain can be swapped.

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
	

	def i(self, ports_v, time=0, get_status=False):
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
		idrain, status = self._get_id(vgs=v1, vgd=v2)
		idrain_np_mos = reversed_factor*sign*idrain

		if self._debug:
                        print "M"+self.descr+":", "vgs:", str(ports_v[0]), "vgd:", str(ports_v[1]), "vds:", str(ports_v[0]-ports_v[1]), status, "id =", str(idrain_np_mos)

		if get_status:
			return idrain_np_mos, status
		else:
			return idrain_np_mos

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
                status = ""
		if vgs <= self.vt and vgd <= self.vt:
			# no channel at both sides
			idr = 0
			status = "off"
		elif vgs > self.vt and vgd > self.vt:
			# zona ohmica: channel at both sides
			idr = .5 * self.kp * (vgs - vgd)*(-2*self.vt + vgs + vgd) * (self.w / self.l)
			status = 'ohmic'
		elif vgs > self.vt and vgd <= self.vt:
			# zona di saturazione: canale al s
			idr =  0.5 * self.kp * ((vgs - self.vt)**2) * (self.w / self.l) * (1 - self.lambd*(vgd - self.vt))
			status = "sat"
		else:
			# zona di saturazione: canale al d
			idr =  -0.5 * self.kp * ((vgd - self.vt)**2) * (self.w / self.l) * (1 - self.lambd*(vgs - self.vt))
			status = "sat"

		return idr, status
	
	def print_op_info(self, ports_v):
		idrain, status = self.i(ports_v, get_status=True)
		print "M"+self.descr+":", status, "vgs:", str(ports_v[0]), "vgd:", str(ports_v[1]), "vds:", str(ports_v[0]-ports_v[1])
		if status == 'sat':
			print "  ", "id =", str(idrain), "gm:", str(self.g(ports_v, 0)), "ro:", str(-1/self.g(ports_v, 1))
		elif status == 'ohmic':
			print "  ", "id =", str(idrain), "gm:", str(self.g(ports_v, 0)), "go:", str(self.g(ports_v, 1))
		else:
			print "  ", "id =", str(idrain), "gm:", str(self.g(ports_v, 0)), "go:", str(self.g(ports_v, 1))
                

if __name__ == '__main__': 
	import numpy
	for type in ('n', 'p'):
		mymos = mosq(nd=2, ng=1, ns=0, kp=1e-3, w=1, l=1, vt=1, lambd=0.5, mos_type=type)
		mymos.descr = "TEST"
		vgs_v = (numpy.array(range(400))-200)/100.0
		vgd_v = [-2, -1.5, -1, -0.5, -0.2, 0, 0.2, 0.5, 1, 1.5, 2]
		for vgd in vgd_v:
			fp = open(type+"mos_test-vgs_sweep-vgd_"+str(vgd)+".txt", 'w')
			fp.write("#VGS\tVGD\tI\tgm\tgos\n")
			for vgs in vgs_v:
				fp.write(str(vgs)+"\t"+str(vgd)+"\t"+str(mymos.i([vgs, vgd]))+"\t"+str(mymos.g([vgs, vgd], 0))+"\t"+str(mymos.g([vgs, vgd], 1))+"\n")
			fp.close()
	
		vgd_v = (numpy.array(range(400))-200)/100.0
		vgs_v = [-2, -1.5, -1, -0.5, -0.2, 0, 0.2, 0.5, 1, 1.5, 2]
		for vgs in vgs_v:
			fp = open(type+"mos_test-vgd_sweep-vgs_"+str(vgs)+".txt", 'w')
			fp.write("#VGD\tVGS\tI\tgos\tgm\n")
			for vgd in vgd_v:
				fp.write(str(vgd)+"\t"+str(vgs)+"\t"+str(mymos.i([vgs, vgd]))+"\t"+str(mymos.g([vgs, vgd], 1))+"\t"+str(mymos.g([vgs, vgd], 0))+"\n")
			fp.close()

	mymos = mosq(nd=2, ng=1, ns=0, kp=1e-3, w=1, l=1, vt=1, lambd=0.5, mos_type='n')
	mymos.descr = "TEST"
	vds_v = range(400)
	vds_v.reverse()
	vds_v = (numpy.array(vds_v))/100.0
	vgs_v = [0.5, 1, 1.2, 1.4, 1.6, 1.8, 2.0]
	fp = open("nmos_test-char.txt", 'w')
	fp.write("#VDS\t")
	mymos.lambd = 0.01
	for vgs in vgs_v:
		fp.write("I(VGS="+str(vgs)+")\t")
	fp.write("\n")
	for vds in vds_v:
		fp.write(str(vds)+"\t")
		for vgs in vgs_v:
			fp.write(str(mymos.i([vgs, vgs-vds]))+"\t")
		fp.write("\n")
	fp.close()
		
	mymos = mosq(nd=2, ng=1, ns=0, kp=1e-3, w=1, l=1, vt=1, lambd=0.5, mos_type='p')
	mymos.descr = "TEST"
	vsd_v = range(400)
	vsd_v.reverse()
	vsd_v = (numpy.array(vsd_v))/100.0
	vgs_v = [-0.5, -1, -1.2, -1.4, -1.6, -1.8, -2.0]
	fp = open("pmos_test-char.txt", 'w')
	fp.write("#VDS\t")
	mymos.lambd = 0.01
	for vgs in vgs_v:
		fp.write("I(VGS="+str(vgs)+")\t")
	fp.write("\n")
	for vsd in vsd_v:
		fp.write(str(vsd)+"\t")
		for vgs in vgs_v:
			fp.write(str(-1*mymos.i([vgs, vgs + vsd]))+"\t")
		fp.write("\n")
	fp.close()
		
