# -*- coding: iso-8859-1 -*-
# switch.py
# Implementation of the voltage controlled switch model
# Copyright 2013 Giuseppe Venturini
# 
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
This module defines two classes:
 switch_device
 switch_model

sn1 o--+         +--o n1
       |         |
      +-+      \ o  
      |R|       \ 
      +-+        +
       |         |
sn2 o--+         +--o n2


"""

import ahkab
import constants, options, utilities, printing
import math 

class switch_device:
	"""This is a general switch element.

	It has the following structure:

	sn1 o--+         +--o n1
	       |         |
	      +-+      \ o  
	      |R|       \ 
	      +-+        +
	       |         |
	sn2 o--+         +--o n2

	The behavior is set by the model supplied.

	The device instance calls the following methods in the model:
	get_i(ports_v, device) - output current
	get_go(ports_v, device) - ouput conductance
	get_gm(ports_v, device) - output transconductance
	get_dc_guess(self, is_on) - guesses for OP

	The device instance accesses the following attributes:
	name, string, the device label
	"""
	def __init__(self, n1, n2, sn1, sn2, model, ic=None):
		"""Voltage controlled switch device
		Parameters:
			n1: output node (+)
			n2: output node (-)
			sn1: input node (+)
			sn2: input node (-)
			model: pass an instance of (v)switch_model
			ic (optional): the initial conditions: True->on, False->off.
		
		Selected methods:
		- get_output_ports() -> (n1, n2)
		- get_drive_ports() -> (n1, n2), (ns1, ns2)
		  

		"""
		class dev_class: pass
		self.device = dev_class()
		self.device.is_on = ic if ic is not None else False
		self.sn1 = sn1
		self.sn2 = sn2
		self.n1 = n1
		self.n2 = n2
		self.ports = ((self.n1, self.n2), (self.sn1, self.sn2))
		self.model = model
		self.opdict = {}
		self.opdict.update({'state':(float('nan'), float('nan'))})
		self.letter_id = 'S'
		self.is_nonlinear = True
		self.is_symbolic = True
		self.dc_guess = self.model.get_dc_guess(self.device.is_on)

	
	def get_drive_ports(self, op):
		"""Get the ports that drive the output ports.

		Returns a tuple of tuples of ports nodes, as:
		(port0, port1, port2...)
		Where each port is in the form:
		port0 = (nplus, nminus)
		"""
		return self.ports

	def get_output_ports(self):
		"""Get the output port, which is (n1,n2) in this case.
		Returns a tuple of tuples of ports nodes, as:
		(port0, port1, port2...)
		Where each port is in the form:
		port0 = (nplus, nminus)
		"""
		return ((self.n1, self.n2),)
	
	def __str__(self):
		rep = self.model.name + " " + str(self.device.is_on)
		return rep

	def i(self, op_index, ports_v, time=0):
		"""Returns the current flowing in the element with the voltages 
		applied	as specified in the ports_v vector.
		
		ports_v: [voltage_across_port0, voltage_across_port1, ...]
		time: the simulation time at which the evaluation is performed. 
		      It has no effect here. Set it to None during DC analysis.
		
		"""
		ret = self.model.get_i(ports_v, self.device)
		# This may be used for debugging
		#print str(ports_v)+" Isw: %g\tRo: %g\tgm: %g" % (ret, 1/self.g(0, ports_v, 0), self.g(0, ports_v, 1))
		return ret

	def update_status_dictionary(self, ports_v):
		"""Updates an internal dictionary that can then be used to provide
		information to the user regarding the status of the element.
		
		Normally, one would call either:
		
		get_op_info(ports_v)
		print_op_info(ports_v)
		"""
		if self.opdict is None:
			self.opdict = {}
		if not (self.opdict['state'] == ports_v[0] and self.opdict.has_key('R')): 
			self.opdict['state'] = ports_v[0]		
			self.opdict['R'] = float(1.0/self.g(0, ports_v[0], 0))
			self.opdict['I'] = float(self.i(0, ports_v[0]))
			self.opdict['STATUS'] = self.device.is_on
		
	def print_op_info(self, ports_v):
		"""Prints out the information regarding the OP status.
		"""
		arr = self.get_op_info(ports_v)
		print arr,

	def get_op_info(self, ports_v):
		"""Operating point info, for design/verification. """
		self.update_status_dictionary(ports_v)

		arr = [["S"+self.descr, 'STATUS:', \
                      "ON"*self.opdict['STATUS']+"OFF"*(not self.opdict['STATUS']), \
                      "VO [V]:", float(self.opdict['state'][0]), \
                      "VS [V]:", float(self.opdict['state'][1]), \
                      "R [ohm]:", self.opdict["R"], \
                      "I [A]:", self.opdict['I'], "",""],]
		#arr.append([  "", "", "", "", "", ""])

		return printing.table_setup(arr)

	
	def g(self, op_index, ports_v, port_index, time=0):
		"""Returns the differential (trans)conductance rs the port specified by port_index
		when the element has the voltages specified in ports_v across its ports,
		at (simulation) time.
		
		ports_v: a list in the form: [voltage_across_port0, voltage_across_port1, ...]
		port_index: an integer, 0 <= port_index < len(self.get_ports())
		time: the simulation time at which the evaluation is performed. Set it to
		None during DC analysis.
		"""
		
		assert op_index == 0 
		assert port_index < 2
		
		if port_index == 0:
			return self.model.get_go(ports_v, self.device)
		if port_index == 1:
			return self.model.get_gm(ports_v, self.device)
		else:
			raise Exception, "Unknown port index passed to switch: bug"

	def get_value_function(self, identifier):
		def get_value(self):
			return self.opdict[identifier]
		return get_value

VT_DEFAULT = 0.0
VH_DEFAULT = 0.0
RON_DEFAULT = 1.
ROFF_DEFAULT = 1./options.gmin

class vswitch_model:
	"""Voltage-controlled switch model.

	sn1 o--+         +--o n1
	       |         |
	      +-+      \ o  
	      |R|       \ 
	      +-+        +
	       |         |
	sn2 o--+         +--o n2


	* R is infinite.
	* The voltage needed to close the switch is: V(sn1)-V(sn2) > VT+VH
	* To re-open it, one needs: V(sn1)-V(sn2) < VT-VH

	The switch commutes between two statuses:

	* ROUT = ROFF
	* ROUT = RON

	None of which can be set to zero or infinite.

	The switching characteristics are modeled with tanh(x).
	"""
	def __init__(self, name, VT=None, VH=None, VON=None, VOFF=None, RON=None, ROFF=None):
		self.name = name
		# convert to VT and VH
		if VON is not None or VOFF is not None:
			VT, VH = self._get_VTVH_from_VONVOFF(float(VON), float(VOFF))
		self.VT = float(VT) if VT is not None else VT_DEFAULT
		self.VH = float(VH) if VH is not None else VH_DEFAULT
		self.RON = float(RON) if RON is not None else RON_DEFAULT
		self.ROFF = float(ROFF) if ROFF is not None else ROFF_DEFAULT
		self.A = (self.RON - self.ROFF)/2
		self.B = (self.RON + self.ROFF)/2.
		self.is_on = False
		self._set_status(self.is_on)
		self.SLOPE = 1e2

	def _get_VTVH_from_VONVOFF(self, VON, VOFF):
		if VON is None or VOFF is None:
			raise CircuitError
		VT = (VON-VOFF)/2.0 + VOFF
		return VT, VT*1e-3

	def _get_V(self, is_on):
		"""Get the effective switching voltage (hyst taken into account)
		"""
		return self.VT + self.VH*2*(not is_on) - self.VH

	def _set_status(self, is_on): 
		"""Set the switch status, which meeans setting the effective 
		switching voltage self.V (w hyst taken into account) 
		"""
		self.V = self.VT + self.VH*2*(not is_on) - self.VH

	def _update_status(self, vin, dev, debug=False):
		"""Check the switch status and move to the other if needed. 
		"""
		Vtest = self._get_V(dev.is_on) 
		R1 = self.A*math.tanh((vin - Vtest)*self.SLOPE) + self.B
		Vtest = self._get_V(not dev.is_on) 
		R2 = self.A*math.tanh((vin - Vtest)*self.SLOPE) + self.B
		self._set_status(dev.is_on) 
		if vin > self.V and not dev.is_on and R1-R2 == 0.0:
			if debug: print "Switching ON: %g" %(vin,)
			dev.is_on = True
			self._set_status(dev.is_on) 
		if vin < self.V and dev.is_on and R1-R2 == 0.0:
			if debug: print "Switching OFF: %g" %(vin,)
			dev.is_on = False
			self._set_status(dev.is_on) 
		self.is_on = dev.is_on
		
	def get_dc_guess(self, is_on):
		"""Returns a list of two floats to be used as initial guesses for the OP analysis
		"""
		return [self.VT*(.9 + is_on*.2)]*2 

	def print_model(self):
		"""All the internal parameters of the model get printed out, 
		for visual inspection. 
		"""
		arr = []
		arr.append([self.name, "", "", "SWITCH MODEL", "", "", "", "", "",  "", "", ""])
		arr.append(["VT", "[V]", self.VT, "VH", "[V]:", self.VH, "RON", "[ohm]", self.RON, "ROFF", "[ohm]", self.ROFF])
		printing.table_print(arr)

	def get_i(self, (vout, vin), dev, debug=False):
		"""Returns the output current.
		"""
		self._update_status(vin, dev)
		R = self.A*math.tanh((vin - self.V)*self.SLOPE) + self.B
		return vout/R
		
	def get_go(self, (vout, vin), dev, debug=False):
		"""Returns the output conductance d(I)/d(Vn1-Vn2)."""
		self._update_status(vin, dev)
		R = self.A*math.tanh((vin - self.V)*self.SLOPE) + self.B
		return 1./R

	def get_gm(self, (vout, vin), dev, debug=False):
		"""Returns the source to output transconductance or d(I)/d(Vsn1-Vsn2)."""
		self._update_status(vin, dev)
		gm = self.A*self.SLOPE*(math.tanh(self.SLOPE*(self.V - vin))**2 - 1)/(self.A*math.tanh(self.SLOPE*(self.V - vin)) - self.B)**2
		return gm+options.gmin

if __name__ == '__main__':
	# This is a small test to check all is OK
	import pylab, numpy, numpy.random
	VT = 10*numpy.random.randn()
	VH = abs(2*numpy.random.randn())
	RON = abs(1e3*numpy.random.randn())
	ROFF = abs(1e4*numpy.random.randn())
	#VT = 0.; VH=1.; RON=100;
	print "Testing a switch with:"
	print "VT: %g\tVH: %g\tRON:%g\tROFF:%g"%(VT, VH, RON, ROFF)
	m = vswitch_model(name='test', VT=VT, VH=VH, RON=RON, ROFF=ROFF)
	VOs = [-12.5, -7.5, -2.5, 2.5, 7.5, 12.5]
	VIs = [-12.5, -7.5, -2.5, 2.5, 7.5, 12.5]
	VMAX = 12.5
	POINTS = 200
	class dev_class: pass
	device = dev_class()
	#TEST 1, VO
	vsweep = (2*VMAX*numpy.arange(POINTS)/float(POINTS)-VMAX)
	vsweep = numpy.concatenate((vsweep, vsweep[::-1]))
	vsweep = vsweep.tolist()
	device.is_on = False
	pylab.hold(True)
	for vin in VIs:
		i = []
		go = []
		gon = []
		for vo in vsweep: 
			i += [m.get_i((vo, vin), device)]
			go += [m.get_go((vo, vin), device)]
			gon += [i[-1]/(vo+1e-12)]
		pylab.subplot(221)
		pylab.plot(vsweep, i, 'o-', label='VIN=%g'%(vin,))
		pylab.ylabel('Output current [A]')
		pylab.xlabel('Output voltage [V]')
		pylab.subplot(222)
		pylab.plot(vsweep, go, 'o', color='g')
		pylab.plot(vsweep, gon, '-')
		pylab.ylabel('Output g [1/V]')
		pylab.xlabel('Output voltage [V]')
	for vo in VOs:
		i = []
		go = []
		gon = []
		for vin in vsweep: 
			i += [m.get_i((vo, vin), device)]
			go += [m.get_go((vo, vin), device)]
			gon += [i[-1]/(vo+1e-12)]
		pylab.subplot(223)
		pylab.plot(vsweep, i, 'o-', label='VO=%g'%(vo,))
		pylab.xlabel('Input voltage [V]')
		pylab.ylabel('Output current [A]')
		pylab.subplot(224)
		pylab.plot(vsweep, go, 'o', color='g', label='VO=%g'%(vo,))
		pylab.plot(vsweep, gon, '-', label='VO=%g'%(vo,))
		pylab.ylabel('Output g [1/V]')
		pylab.xlabel('Input voltage [V]')
	pylab.subplot(221)
	pylab.legend()
	pylab.subplot(222)
	pylab.legend()
	pylab.subplot(223)
	pylab.legend()
	pylab.subplot(224)
	pylab.legend()
	pylab.show()
