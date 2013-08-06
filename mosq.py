# -*- coding: iso-8859-1 -*-
# mosq.py
# Implementation of the square-law MOS transistor model
# Copyright 2012 Giuseppe Venturini
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
 mosq_device
 mosq_model

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

"""

import constants, options, utilities, printing
import math 


# DEFAULT VALUES FOR 500n CH LENGTH
COX_DEFAULT = .7e-3
VTO_DEFAULT = .5
GAMMA_DEFAULT = 1
PHI_DEFAULT = .7
KP_DEFAULT = 50e-6
LAMBDA_DEFAULT = .5
AVT_DEFAULT = 7.1e-3*1e-6
AKP_DEFAULT = 1.8e-2*1e-6

TCV_DEFAULT = 1e-3
BEX_DEFAULT = -1.5 

ISMALL_GUESS_MIN = 1e-10

class mosq_device:
	def __init__(self, nd, ng, ns, nb, W, L, model, M=1, N=1):
		""" EKV device
		Parameters:
			nd: drain node
			ng: gate node
			ns: source node
			nb: bulk node
			L: element width [m]
			W: element length [m]
			M: multiplier (n. of shunt devices)
			N: series mult. (n. of series devices)
			model: pass an instance of ekv_mos_model
		
		Selected methods:
		- get_output_ports() -> (nd, ns)
		- get_drive_ports() -> (nd, nb), (ng, nb), (ns, nb)
		  

		"""
		self.ng = ng
		self.nb = nb
		self.n1 = nd
		self.n2 = ns
		self.ports = ((self.n1, self.n2), (self.ng, self.n2), (self.nb, self.n2))
		class dev_class: pass # empty class to hold device parameters
		self.device = dev_class()
		self.device.L = float(L) #channel length -
		self.device.W = float(W) #channel width -
		self.device.M = int(M) #parallel multiple device number
		self.device.N = int(N) #series multiple device number
		self.device.mckey = None
		self.mosq_model = model
		self.mc_enabled = False
		self.opdict = {}
		self.opdict.update({'state':(float('nan'), float('nan'), float('nan'))})
		self.letter_id = 'M'
		self.is_nonlinear = True
		self.dc_guess = [self.mosq_model.VTO*(0.4)*self.mosq_model.NPMOS, self.mosq_model.VTO*(1.1)*self.mosq_model.NPMOS, 0]

		devcheck, reason =  self.mosq_model._device_check(self.device)
		if not devcheck:
			raise Exception, reason + " out of boundaries."		
	
	def get_drive_ports(self, op):
		"""Returns a tuple of tuples of ports nodes, as:
		(port0, port1, port2...)
		Where each port is in the form:
		port0 = (nplus, nminus)
		"""
		return self.ports #d,g,b

	def get_output_ports(self):
		return ((self.n1, self.n2),)
	
	def __str__(self):
		mos_type = self._get_mos_type()
		rep = " " + self.mosq_model.name + " w="+ str(self.device.W) + " l=" + \
		str(self.device.L) + " M="+ str(self.device.M) + " N=" + \
		str(self.device.N)

		return rep

	def _get_mos_type(self):
		"""Returns N or P (capitalized)
		"""
		mtype = 'N' if self.mosq_model.NPMOS == 1 else 'P'
		return mtype

	def i(self, op_index, ports_v, time=0):
		"""Returns the current flowing in the element with the voltages 
		applied	as specified in the ports_v vector.
		
		ports_v: [voltage_across_port0, voltage_across_port1, ...]
		time: the simulation time at which the evaluation is performed. 
		      It has no effect here. Set it to None during DC analysis.
		
		"""
		#print ports_v
		ret = self.mosq_model.get_ids(self.device, ports_v, self.opdict)
		
		return ret

	def update_status_dictionary(self, ports_v):
		if self.opdict is None:
			self.opdict = {}
		if not (self.opdict['state'] == ports_v[0] and self.opdict.has_key('gmd')) or \
			not (self.opdict['state'] == ports_v[0] and self.opdict.has_key('gm')) or \
			not (self.opdict['state'] == ports_v[0] and self.opdict.has_key('gmb')) or \
			not (self.opdict['state'] == ports_v[0] and self.opdict.has_key('Ids')):

			self.opdict['state'] == ports_v[0]		
			self.opdict['gmd'] = self.g(0, ports_v[0], 0)
			self.opdict['gm'] = self.g(0, ports_v[0], 1)
			self.opdict['gmb'] = self.g(0, ports_v[0], 2)
			self.opdict['Ids'] = self.i(0, ports_v[0])
		
	def print_op_info(self, ports_v):
		arr = self.get_op_info(ports_v)
		print arr,

	def get_op_info(self, ports_v):
		"""Operating point info, for design/verification. """
		mos_type = self._get_mos_type()
		self.update_status_dictionary(ports_v)
		sat_status = "SATURATION" if self.opdict['SAT'] else "LINEAR"
		if not self.opdict["ON"]: 
			status = "OFF"
		else:
			status = "ON"

		arr = [["M"+self.descr, mos_type.upper()+" ch", status, "", "", sat_status, "", "", "", "", "",""],]
		arr.append(["beta", "[A/V^2]:", self.opdict['beta'], "Weff", "[m]:", str(self.opdict['W'])+" ("+str(self.device.W)+")", "L", "[m]:", str(self.opdict['L'])+ " ("+str(self.device.L)+")", "M/N:", "", str(self.device.M)+"/"+str(self.device.N)])
		arr.append(["Vds", "[V]:", float(ports_v[0][0]), "Vgs", "[V]:", float(ports_v[0][1]), "Vbs", "[V]:", float(ports_v[0][2]),  "", "", ""])
		arr.append([ "VTH", "[V]:", self.opdict['VTH'], "VOD", "[V]:", self.opdict['VOD'], "", "","", "VA", "[V]:", str(self.opdict['Ids']/self.opdict['gmd'])])
		arr.append(["Ids", "[A]:", self.opdict['Ids'], "", "", "", "", "", "", "", "", '']) 
		arr.append(["gm", "[S]:", self.opdict['gm'], "gmb", "[S]:", self.opdict['gmb'], "ro", "[Ohm]:", 1/self.opdict['gmd'], "", "", ""])
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
		assert port_index < 3
		
		if port_index == 0:
			g = self.mosq_model.get_gmd(self.device, ports_v, self.opdict)
		elif port_index == 1:
			g = self.mosq_model.get_gm(self.device, ports_v, self.opdict)
		if port_index == 2:
			g = self.mosq_model.get_gmb(self.device, ports_v, self.opdict)

		if op_index == 0 and g == 0:
			if port_index == 2:
				sign = -1
			else:
				sign = +1
			g = sign*options.gmin*2

		#print type(g), g

		if op_index == 0 and port_index == 0:
			self.opdict.update({'gmd':g})
		elif op_index == 0 and port_index == 1:
			self.opdict.update({'gm':g})
		elif op_index == 0 and port_index == 2:
			self.opdict.update({'gmb':g})

		return g

	def get_value_function(self, identifier):
		def get_value(self):
			return self.opdict[identifier]
		return get_value

	def get_mc_requirements(self):
		return True, 2
	def setup_mc(self, status, mckey):
		self.mc_enabled = status
		if self.mc_enabled:
			self.device.mckey = mckey
		else:
			self.device.mckey = None

class scaling_holder: pass # will hold the scaling factors

class mosq_mos_model:
	def __init__(self, name=None, TYPE='n', TNOM=None, COX=None, \
	GAMMA=None, NSUB=None, PHI=None, VTO=None, KP=None, \
	LAMBDA=None, AKP=None, AVT=None,\
	TOX=None, VFB=None, U0=None, TCV=None, BEX=None):
		
		self.scaling = scaling_holder()

		self.name = "model_mosq0" if name is None else name
		Vth = constants.Vth()
		self.TNOM = float(TNOM) if TNOM is not None else constants.Tref	
		#print "TYPE IS:" + TYPE		
		self.NPMOS = 1 if TYPE == 'n' else -1

		# optional parameters (no defaults)	
		self.TOX = float(TOX) if TOX is not None else None
		self.NSUB = float(NSUB)  if NSUB is not None else None
		self.VFB = self.NPMOS*float(VFB) if VFB is not None else None
		self.U0 = float(U0) if U0 is not None else None

		# crucial parameters
		if COX is not None: 
			self.COX = float(COX)
		elif TOX is not None:
			self.COX = constants.si.eox/TOX
		else:
			self.COX = COX_DEFAULT
	
		if GAMMA is not None:
			self.GAMMA = float(GAMMA)
		elif NSUB is not None:
			self.GAMMA = math.sqrt(2*constants.e*constants.si.esi*NSUB*10**6/self.COX)
		else:
			self.GAMMA = GAMMA_DEFAULT
		if PHI is not None:
			self.PHI = float(PHI)
		elif NSUB is not None:
			self.PHI = 2*constants.Vth(self.TNOM)*math.log(NSUB*10**6/constants.si.ni(self.TNOM))
		else:
			self.PHI = PHI_DEFAULT
		if VTO is not None:
			self.VTO = self.NPMOS*float(VTO)
			if self.VTO < 0:
				print "(W): model %s has internal negative VTO (%f V)." % (self.name, self.VTO)
		elif VFB is not None:
			self.VTO = VFB + PHI + GAMMA*PHI #inv here??
		else:
			self.VTO = self.NPMOS*VTO_DEFAULT

		if KP is not None:
			self.KP = float(KP)
		elif U0 is not None:
			self.KP = (U0*10**-4)*self.COX
		else:
			self.KP = KP_DEFAULT

		self.LAMBDA = LAMBDA if LAMBDA is not None else LAMBDA_DEFAULT
		# Intrinsic model temperature parameters
		self.TCV = self.NPMOS*float(TCV) if TCV is not None else self.NPMOS*TCV_DEFAULT
		self.BEX = float(BEX) if BEX is not None else BEX_DEFAULT
	
		# Monte carlo
		self.AVT = AVT if AVT is not None else AVT_DEFAULT
		self.AKP = AKP if AKP is not None else AKP_DEFAULT

		self.set_device_temperature(constants.T)

		sc, sc_reason = self._self_check()
		if not sc:
			raise Exception, sc_reason + " out of range"

	def set_device_temperature(self, T):
		"""Change the temperature of the device. VTO, KP and PHI get updated.
		"""		
		self.TEMP = T
		self.VTO = self.VTO - self.TCV*(T-self.TNOM)
		self.KP = self.KP*(T/self.TNOM)**self.BEX
		self.PHI = self.PHI * T/self.TNOM + 3*constants.Vth(self.TNOM)*math.log(T/self.TNOM) \
			   - constants.si.Eg(self.TNOM)*T/self.TNOM + constants.si.Eg(T)

	def get_device_temperature(self):
		"""Returns the temperature of the device - in K.
		"""
		return self.TEMP

	def print_model(self):
		"""All the internal parameters of the model get printed out, 
		for visual inspection. Notice some can be set to None 
		(ie not available) if they were not provided in the netlist
		or some not provided are calculated from the others.
		"""
		arr = []
		TYPE = 'N' if self.NPMOS == 1 else "P"
		arr.append([self.name, "", "", TYPE+" MOS", "SQUARE MODEL", "", "", "", "",  "", "", ""])
		arr.append(["KP", "[A/V^2]", self.KP, "VTO", "[V]:", self.VTO, "TOX", "[m]", self.TOX, "COX", "[F/m^2]:", self.COX])
		arr.append(["PHI", "[V]:", self.PHI, "GAMMA", "sqrt(V)", self.GAMMA, "NSUB", "[cm^-3]", self.NSUB,  "VFB", "[V]:", self.VFB])
		arr.append(["U0", "[cm^2/(V*s)]:", self.U0, "TCV", "[V/K]", self.TCV, "BEX", "", self.BEX,  "", "", ""])
		printing.table_print(arr)

	def get_voltages(self, vds, vgs, vbs):
		"""Performs the D <-> S swap if needed.
		Returns: 
		(VDS, VGS, VBS) after the swap
		CS, an integer which equals to:
  			+1 if no swap was necessary,
			-1 if VD and VS have been swapped.
		"""
		# vd / vs swap
		vds = float(vds)
		vgs = float(vgs)
		vbs = float(vbs)
		vds = vds*self.NPMOS
		vgs = vgs*self.NPMOS
		vbs = vbs*self.NPMOS
		if vds < 0:
			vds_new = -vds
			vgs_new = vgs - vds
			vbs_new = vbs - vds
			cs = -1
		else:
			vds_new = vds
			vgs_new = vgs
			vbs_new = vbs
			cs = +1
		#print ((float(vds_new), float(vgs_new), float(vbs_new)), cs)
		return ((float(vds_new), float(vgs_new), float(vbs_new)), cs)

	def get_svt_skp(self, device, debug=False):
		if device.mckey and debug:
			print "Monte carlo enabled. key:", device.mckey
		if device.mckey:
			svt = device.mckey[0] * self.AVT / math.sqrt(2*device.W*device.L)
			skp = device.mckey[1] * self.AKP / math.sqrt(2*device.W*device.L)
		else:
			svt, skp = 0, 0
		return svt, skp

	def get_ids(self, device, (vds, vgs, vbs), opdict=None, debug=False):
		"""Returns:
			IDS, the drain-to-source current (de-normalized),
			qs, the (scaled) charge at the source,
			qr, the (scaled) charge at the drain.
		"""
		if debug: 
			print "=== Current for vds:", vds, "vgs:", vgs, "vbs:", vbs
		
		(vds, vgs, vbs), CS_FACTOR = self.get_voltages(vds, vgs, vbs)

		# monte carlo support		
		svt, skp = self.get_svt_skp(device, debug=debug)

		#print "PHI:", self.PHI, "vbs:", vbs
		VT = self.VTO + svt + self.GAMMA*(math.sqrt(-vbs+2*self.PHI) - math.sqrt(2*self.PHI))
		if vgs < VT:
			ids = options.iea*(vgs/VT+vds/VT)/100 
		else:
			if vds < vgs - VT:
				ids = (skp+1)*self.KP*device.W/device.L*((vgs-VT)*vds - .5*vds**2)
			else:
				ids = (skp+1)*.5*self.KP*device.W/device.L*(vgs-VT)**2*(1+self.LAMBDA*(vds-vgs+VT))

		Ids =  CS_FACTOR * self.NPMOS* device.M/device.N * ids

		vds_real = vds if CS_FACTOR == 1 else -vds
		vgs_real = vgs if CS_FACTOR == 1 else vgs-vds
		vbs_real = vbs if CS_FACTOR == 1 else vbs-vds

		opdict.update({'state':(vds_real*self.NPMOS, vgs_real*self.NPMOS, vbs_real*self.NPMOS)})
		opdict.update({'Ids':Ids, "W":device.W, "L":device.L, "ON":1*(vgs>=VT)})
		opdict.update({'beta':.5*self.KP*device.W/device.L})
		opdict.update({'VTH':VT, "VOD":self.NPMOS*(vgs-VT), 'SAT':vds>vgs-VT})


		return Ids 

		#if debug:
		#	print "vd:", vd, "vg:",VG/self.scaling.Ut, "vs:", vs, "vds:", vd-vs
		#	print "v_ifn:", v_ifn, "v_irn:",v_irn
		#	print "ifn:", ifn, "irn:",irn
		#	print "ip_abs_err:", ip_abs_err
		#	print "Vth:", self.scaling.Ut
		#	print "nv", nv, "Is", self.scaling.Is
		#	print "Weff:", device.W, "Leff:", Leff
		#	print "NPMOS:", self.NPMOS, "CS_FACTOR", CS_FACTOR


		
	def get_gmb(self, device, (vds, vgs, vbs), opdict=None, debug=False):
		"""Returns the source-bulk transconductance or d(IDS)/d(VS-VB)."""
		svt, skp = self.get_svt_skp(device, debug=False)
		(vds, vgs, vbs), CS_FACTOR = self.get_voltages(vds, vgs, vbs)
		VT = self.VTO + svt + self.GAMMA*(math.sqrt(-vbs+2*self.PHI) - math.sqrt(2*self.PHI))
		if CS_FACTOR < 0:
			return CS_FACTOR*self.NPMOS*self.get_gmb(device, (vds*self.NPMOS, vgs*self.NPMOS, vbs*self.NPMOS), opdict, debug)
		if vgs < VT:
			gmb = 0 
		else:
			if vds < vgs - VT:
				gmb = self.KP*self.GAMMA*vds*device.W/(2*device.L*(2*self.PHI - vbs)**(1.0/2))
			else:
				gmb = -0.25*self.KP*self.GAMMA*self.LAMBDA*device.W*\
				(-self.GAMMA*(-2**(1.0/2)*self.PHI**(1.0/2) + (2*self.PHI - vbs)**(1.0/2)) + vgs - self.VTO)**2/(device.L*(2*self.PHI - vbs)**(1.0/2)) \
				+ 0.5*self.KP*self.GAMMA*device.W*(self.LAMBDA*(self.GAMMA*(-2**(1.0/2)*self.PHI**(1.0/2) + (2*self.PHI - vbs)**(1.0/2)) + vds - vgs + self.VTO) + 1.0)*\
				(-self.GAMMA*(-2**(1.0/2)*self.PHI**(1.0/2) + (2*self.PHI - vbs)**(1.0/2)) + vgs - self.VTO)/(device.L*(2*self.PHI - vbs)**(1.0/2))
		gmb = self.NPMOS * (1+skp) * gmb * device.M/device.N
		return gmb

	def get_gmd(self, device, (vds, vgs, vbs), opdict=None, debug=False):
		"""Returns the drain-bulk transconductance or d(IDS)/d(VD-VB)."""
		svt, skp = self.get_svt_skp(device, debug=False)
		(vds, vgs, vbs), CS_FACTOR = self.get_voltages(vds, vgs, vbs)
		VT = self.VTO + svt + self.GAMMA*(math.sqrt(-vbs+2*self.PHI) - math.sqrt(2*self.PHI))
		if vgs < VT:
			gmd = options.iea/VT/100 
		else:
			if vds < vgs - VT:
				gmd = self.KP*device.W/device.L*(-self.GAMMA*(-2**(1.0/2)*self.PHI**(1.0/2) + (2*self.PHI - vbs)**(1.0/2)) - 1.0*vds + vgs - self.VTO)
			else:
				gmd = 0.5*self.KP*self.LAMBDA*device.W/device.L*(-self.GAMMA*(-2**(1.0/2)*self.PHI**(1.0/2) + (2*self.PHI - vbs)**(1.0/2)) + vgs - self.VTO)**2
		gmd = (1+skp) * gmd * device.M/device.N
		return gmd

	def get_gm(self, device, (vds, vgs, vbs), opdict=None, debug=False):
		"""Returns the gate-bulk transconductance or d(IDS)/d(VG-VB)."""
		svt, skp = self.get_svt_skp(device, debug=False)
		(vds, vgs, vbs), CS_FACTOR = self.get_voltages(vds, vgs, vbs)
		if CS_FACTOR < 0:
			return self.get_gm(device, (vds*self.NPMOS, vgs*self.NPMOS, vbs*self.NPMOS), opdict, debug)
		VT = self.VTO + svt + self.GAMMA*(math.sqrt(-vbs+2*self.PHI) - math.sqrt(2*self.PHI))
		if vgs < VT:
			gm = options.iea/VT/100 
		else:
			if vds < vgs - VT:
				gm =  self.KP*device.W/device.L*vds
			else:
				gm = -0.5*self.KP*self.LAMBDA*device.W/device.L*(-self.GAMMA*(-2**(1.0/2)*self.PHI**(1.0/2) + (2*self.PHI - vbs)**(1.0/2)) + vgs - self.VTO)**2 \
					 + 0.5*self.KP*device.W/device.L*(self.LAMBDA*(self.GAMMA*(-2**(1.0/2)*self.PHI**(1.0/2) + (2*self.PHI - vbs)**(1.0/2)) + vds - vgs + self.VTO) + 1.0)*\
					(-2*self.GAMMA*(-2**(1.0/2)*self.PHI**(1.0/2) + (2*self.PHI - vbs)**(1.0/2)) + 2*vgs - 2*self.VTO)
		gm = CS_FACTOR * self.NPMOS * (1+skp) * gm * device.M/device.N
		return gm

	def _self_check(self):
		"""Performs sanity check on the model parameters."""
		ret = True, ""
		if self.NSUB is not None and self.NSUB < 0:
			ret = (False, "NSUB "+str(self.NSUB))
		elif self.U0 is not None and not self.U0 > 0:
			ret = (False, "UO "+str(self.U0))
		elif not self.GAMMA > 0:			 
			ret = (False, "GAMMA "+str(self.GAMMA))
		elif not self.PHI > 0.1:			 
			ret = (False, "PHI "+str(self.PHI))
		elif self.AVT and self.AVT < 0:
			ret = (False, "AVT "+str(self.AVT))
		elif self.AKP and self.AKP < 0:
			ret = (False, "AKP "+str(self.AKP))
		return ret

	def _device_check(self, adev):
		"""Performs sanity check on the device parameters."""
		if not adev.L > 0:
			ret = (False, "L")
		elif not adev.W > 0:
			ret = (False, "W")
		elif not adev.N > 0:
			ret = (False, "N")
		elif not adev.M > 0:
			ret = (False, "M")
		else:
			ret = (True, "")
		return ret
		
if __name__ == '__main__':
	# Tests
	import matplotlib.pyplot as plt
	import numpy

	m = mosq_mos_model(TYPE='p', KP=50e-6, VTO=.4)
	ma = mosq_device(1, 2, 3, 4, W=10e-6,L=1e-6, model=m)
	ma.descr = "1"

	# OP test
	vds = numpy.arange(0, 100)/100.0*5-2.5
	vgs = -.55
	vbs = 2
	#ma.print_op_info(((vds, vgs, vbs),))
	#m.print_model()
	i= []
	g=[]
	for X in vds:
		i+=[ma.i(0, (X,  vgs, vbs))]
		g += [ma.g(0, (X, vgs, vbs), 0)]
	plt.figure()
	plt.plot(vds, i)
	plt.hold(True)
	plt.plot(vds, g)
	gart = (numpy.array(i[1:]) - numpy.array(i[:-1]))/(vds[1]-vds[0])
	plt.plot(vds[1:], gart)
	plt.show()
