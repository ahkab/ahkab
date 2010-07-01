# -*- coding: iso-8859-1 -*-
# ekv.py
# Partial implementation of the EKV 3.0 MOS transistor model
# Copyright 2010 Giuseppe Venturini
# 
# The EKV model was developed by Matthias Bucher, Christophe Lallement, 
# Christian Enz, Fabien Théodoloz, François Krummenacher at the Electronics 
# Laboratories, Swiss Federal Institute of Technology (EPFL), 
# Lausanne, Switzerland. 
# This implementation is based upon:
# 1. Matthias Bucher, Christian Enz, François Krummenacher, Jean-M. Sallese, 
# Christophe Lallement and Alain-S. Porret, 
# The EKV 3.0 Compact MOS Transistor Model: Accounting for Deep-Submicron 
# Aspects, <http://www.nsti.org/publications/MSM/2002/pdf/346.pdf>
# 2. EKV 2.6 Technical report, <http://legwww.epfl.ch/ekv/pdf/ekv_v262.pdf>.
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
The EKV model was developed by Matthias Bucher, Christophe Lallement, 
Christian Enz, Fabien Théodoloz, François Krummenacher at the Electronics 
Laboratories, Swiss Federal Institute of Technology (EPFL), 
Lausanne, Switzerland. The Tecnical Report upon which this implementation 
is based is available here:
<http://legwww.epfl.ch/ekv/pdf/ekv_v262.pdf>.

This module defines two classes:
 ekv_device
 ekv_mos_model


Features:
- EKV model implementation, computation of charges, potentials,
  reverse and forward currents, slope factor and normalization factors,
- Calculation of trans-conductances based on the charge approach,
- Rudimentary temperature effects.

The Missing Features:
- Channel length modulation
- Reverse Short Channel Effect (RSCE)


TODO:
- Complex mobility degradation is missing
- Transcapacitances
- Quasistatic implementation
"""

import constants, options, utilities, printing
import math 


# DEFAULT VALUES FOR 500n CH LENGTH
COX_DEFAULT = .7e-3

VTO_DEFAULT = .5
GAMMA_DEFAULT = 1
PHI_DEFAULT = .7
KP_DEFAULT = 50e-6

TCV_DEFAULT = 1e-3
BEX_DEFAULT = -1.5 

class ekv_device:
	INIT_IFRN_GUESS = 1e-3
	def __init__(self, nd, ng, ns, nb, W, L, model, M=1, N=1):
		""" ekv device
		nd: drain node
		ng: gate node
		ns: source node
		nb: bulk node
		L: element width [m]
		W: element length [m]
		M: multiplier (n. of shunt devices)
		N: series mult. (n. of series devices)
		model: pass an instance of ekv_mos_model
		"""
		self.ng = ng
		self.nb = nb
		self.n1 = nd
		self.n2 = ns
		self.ports = ((self.n1, self.nb), (self.ng, self.nb), (self.n2, self.nb))
		class dev_class: pass # empty class to hold device parameters
		self.device = dev_class()
		self.device.L = float(L) #channel length -
		self.device.W = float(W) #channel width -
		self.device.M = int(M) #parallel multiple device number
		self.device.N = int(N) #series multiple device number
		self.ekv_model = model
		self.opdict = {}
		self.opdict.update({'state':(float('nan'), float('nan'), float('nan'))})
		self.opdict.update({'ifn':self.INIT_IFRN_GUESS})
		self.opdict.update({'irn':self.INIT_IFRN_GUESS})
		self.opdict.update({'ip_abs_err':self.ekv_model.get_ip_abs_err(self.device)})
		self.letter_id = 'M'
		self.is_nonlinear = True
		self.dc_guess = [self.ekv_model.VTO*(0.1)*self.ekv_model.NPMOS, self.ekv_model.VTO*(1.1)*self.ekv_model.NPMOS, 0]

		devcheck, reason =  self.ekv_model._device_check(self.device)
		if not devcheck:
			raise Exception, reason + " out of boundaries."		
	
	def get_drive_ports(self, op):
		"""Returns a tuple of tuples of ports nodes, as:
		(port0, port1, port2...)
		Where each port is in the form:
		port0 = (nplus, nminus)
		"""
		return self.ports #d,g,s

	def get_output_ports(self):
		return ((self.n1, self.n2),)
	
	def __str__(self):
		mos_type = self._get_mos_type()
		rep = "type=" + mos_type + " w="+ str(self.device.W) + " l=" + \
		str(self.device.L) + " M="+ str(self.device.M) + " N=" + \
		str(self.device.N)

		return rep

	def _get_mos_type(self):
		"""Returns N or P (capitalized)
		"""
		mtype = 'N' if self.ekv_model.NPMOS == 1 else 'P'
		return mtype

	def i(self, op_index, ports_v, time=0):
		"""Returns the current flowing in the element with the voltages 
		applied	as specified in the ports_v vector.
		
		ports_v: [voltage_across_port0, voltage_across_port1, ...]
		time: the simulation time at which the evaluation is performed. 
		      Set it to None during DC analysis.
		
		"""
		#self.opdict = {'state':ports_v[0]}
		#print ports_v
		ret, j1, j2 = self.ekv_model.get_ids(self.device, ports_v, \
				self.opdict)
		
		return ret

	def update_status_dictionary(self, ports_v):
		if self.opdict is None:
			self.opdict = {}
		if not (self.opdict['state'] == ports_v[0] and self.opdict.has_key('gmd')) or \
			not (self.opdict['state'] == ports_v[0] and self.opdict.has_key('gmg')) or \
			not (self.opdict['state'] == ports_v[0] and self.opdict.has_key('gms')) or \
			not (self.opdict['state'] == ports_v[0] and self.opdict.has_key('Ids')):

			self.opdict['state'] == ports_v[0]		
			self.opdict['gmd'] = self.g(0, ports_v[0], 0)
			self.opdict['gmg'] = self.g(0, ports_v[0], 1)
			self.opdict['gms'] = self.g(0, ports_v[0], 2)
			self.opdict['Ids'] = self.i(0, ports_v[0])

		gmd = self.opdict['gmd']
		gmg = self.opdict['gmg']		
		gms = self.opdict['gms']		
		ids = self.opdict['Ids']

		if ids == 0:
			TEF = float('nan')
		else:
			TEF = -gms*constants.Vth()/ids
		self.opdict['TEF'] = TEF
		

			
	def print_op_info(self, ports_v):
		"""Operating point info, for design/verification. """
		mos_type = self._get_mos_type()

		self.update_status_dictionary(ports_v)

		sat_status = "SATURATION" if self.opdict['SAT'] else "LINEAR"
		if self.opdict["WMSI"] == 0: 
			wmsi_status = "WEAK INVERSION"
		if self.opdict["WMSI"] == 1: 
			wmsi_status = "MODERATE INVERSION"
		if self.opdict["WMSI"] == 2: 
			wmsi_status = "STRONG INVERSION"

		arr = [["M"+self.descr, mos_type.upper()+" ch",wmsi_status, "", "", sat_status, "", "", "", "", "",""],]
		arr.append(["beta", "[A/V^2]:", self.opdict['beta'], "Weff", "[m]:", str(self.opdict['Weff'])+" ("+str(self.device.W)+")", "Leff", "[m]:", str(self.opdict['Leff'])+ " ("+str(self.device.L)+")", "M/N:", "", str(self.device.M)+"/"+str(self.device.N)])
		arr.append(["Vdb", "[V]:", float(ports_v[0][0]), "Vgb", "[V]:", float(ports_v[0][1]), "Vsb", "[V]:", float(ports_v[0][2]),  "Vp", "[V]:", self.opdict['Vp'],])
		arr.append([ "VTH", "[V]:", self.opdict['VTH'], "VOD", "[V]:", self.opdict['VOD'], "", "", "", "VA", "[V]:", str(self.opdict['Ids']/self.opdict['gmd'])])
		arr.append(["Ids", "[A]:", self.opdict['Ids'], "n: ", "",self.opdict['n'], "Ispec", "[A]:", self.opdict["Ispec"], "TEF:", "", str(self.opdict['TEF']),]) 
		arr.append(["gmg", "[S]:", self.opdict['gmg'], "gms", "[S]:", self.opdict['gms'], "rob", "[Ohm]:", 1/self.opdict['gmd'], "", "", ""])
		arr.append(["if:", "", self.opdict['ifn'],"ir:", "", self.opdict['irn'], "Qf", "[C/m^2]:", self.opdict["qf"], "Qr", "[C/m^2]:", self.opdict["qf"],])
		#arr.append([  "", "", "", "", "", ""])

		printing.table_print(arr)

	
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
			g = self.ekv_model.get_gmd(self.device, ports_v, self.opdict)
		elif port_index == 1:
			g = self.ekv_model.get_gmg(self.device, ports_v, self.opdict)
		if port_index == 2:
			g = self.ekv_model.get_gms(self.device, ports_v, self.opdict)

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
			self.opdict.update({'gmg':g})
		elif op_index == 0 and port_index == 2:
			self.opdict.update({'gms':g})

		return g

	def get_value_function(self, identifier):
		def get_value(self):
			return self.opdict[identifier]
		return get_value

class scaling_holder: pass # will hold the scaling factors

class ekv_mos_model:
	def __init__(self, name=None, TYPE='n', TNOM=None, COX=None, \
	GAMMA=None, NSUB=None, PHI=None, VTO=None, KP=None, \
	TOX=None, VFB=None, U0=None, TCV=None, BEX=None):
		
		self.scaling = scaling_holder()

		self.name = "model_ekv0" if name is None else name
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
			self.COX = constants.si.esi/TOX
		else:
			self.COX = COX_DEFAULT
	
		if GAMMA is not None:
			self.GAMMA = float(GAMMA)
		elif NSUB is not None:
			math.sqrt(2*constants.e*constants.si.esi*NSUB*10**6/self.COX)
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
	
		# Intrinsic model temperature parameters
		self.TCV = self.NPMOS*float(TCV) if TCV is not None else self.NPMOS*TCV_DEFAULT
		self.BEX = float(BEX) if BEX is not None else BEX_DEFAULT
	
		self.set_device_temperature(constants.T)

		#Setup switches
		self.SATLIM = math.exp(4)
		self.WMSI_factor = 10
		self.NQS = 0
		self.CLM_SWITCH = True

		sc, sc_reason = self._self_check()
		if not sc:
			raise Exception, sc_reason + " out of range"

	def set_device_temperature(self, T):
		self.TEMP = T
		self.VTO = self.VTO - self.TCV*(T-self.TNOM)
		self.KP = self.KP*(T/self.TNOM)**self.BEX
		self.PHI = self.PHI * T/self.TNOM + 3*constants.Vth(self.TNOM)*math.log(T/self.TNOM) \
			   - constants.si.Eg(self.TNOM)*T/self.TNOM + constants.si.Eg(T)


	def get_device_temperature(self):
		return self.TEMP

	def get_voltages(self, vd, vg, vs):
		# vd / vs swap
		vd = vd*self.NPMOS
		vg = vg*self.NPMOS
		vs = vs*self.NPMOS
		if vs > vd:
			vd_new = vs
			vs_new = vd
			cs = -1
		else:
			vd_new = vd
			vs_new = vs
			cs = +1
		#print "got", vd,vg,vs,"returned", vd_new, vg, vs_new
		return ((float(vd_new), float(vg), float(vs_new)), cs)

	def get_ip_abs_err(self, device):
		return  options.iea / (4*constants.Vth()**2*self.KP*device.M*device.W/device.L)

	def setup_scaling(self, nq, device, debug=False):
		self.scaling.Ut = constants.Vth()
		self.scaling.Is = 2 * nq * self.scaling.Ut**2 * self.KP * device.W/device.L
		self.scaling.Gs = 2 * nq * self.scaling.Ut * self.KP * device.W/device.L
		self.scaling.Qs = 2 * nq * self.scaling.Ut * self.COX
		
		if debug:
			print "Is", self.scaling.Is
		return		
	
	def get_vp_nv_nq(self, VG):
		#VGeff = VG - self.VTO + self.PHI + self.GAMMA*math.sqrt(self.PHI)
		##if VGeff > 0:		
		#	VP = VGeff - self.PHI - self.GAMMA*(math.sqrt(VGeff + (self.GAMMA/2)**2) -self.GAMMA/2)
		#else:
		#	VP = -self.PHI
		VGeff = VG - self.VTO + self.PHI + self.GAMMA*math.sqrt(self.PHI)
		if VGeff > 0:		
			VP = VG - self.VTO - self.GAMMA*(math.sqrt(VG -self.VTO +(math.sqrt(self.PHI)+self.GAMMA/2)**2) -(math.sqrt(self.PHI)+self.GAMMA/2))
		else:
			VP = -self.PHI		
		

		nq = 1 + .5 * self.GAMMA / math.sqrt(self.PHI + .5*VP)
		nv = 1 + .5 * self.GAMMA / math.sqrt(self.PHI +    VP)

		return VP, nv, nq

	def get_ids(self, device, (vd, vg, vs), opdict=None, debug=False):
		if debug: print "Current for vd:", vd, "vg:", vg, "vs:", vs
		ip_abs_err = self.get_ip_abs_err(device) if opdict['ip_abs_err'] is None else opdict['ip_abs_err']
		
		(VD, VG, VS), CS_FACTOR = self.get_voltages(vd, vg, vs)	

		#Weff, Leff = self.get_eff_wl(device.W, device.L)		

		VP, nv, nq = self.get_vp_nv_nq(VG)
		
		self.setup_scaling(nq, device)

		vp = VP/self.scaling.Ut
		vs = VS/self.scaling.Ut
		vd = VD/self.scaling.Ut

		if debug: print "Scaled voltages: vd:", vd, "vp:", vp, "vs:", vs

		v_ifn = vp - vs
		ifn = self.get_ismall(v_ifn, opdict['ip_abs_err'], opdict['ifn'])

		Leff = device.L

		v_irn = vp - vd
		irn = self.get_ismall(v_irn, opdict['ip_abs_err'], opdict['irn'])
		
		if debug:
			print "vd:", vd, "vg:",vg/self.scaling.Ut, "vs:", vs, "vds:", vd-vs
			print "v_ifn:", v_ifn, "v_irn:",v_irn
			print "ip_abs_err:", ip_abs_err
			print "Vth:", self.scaling.Ut
			print "Weff:", device.W, "Leff:", Leff
			print "NPMOS:", self.NPMOS, "CS_FACTOR", CS_FACTOR

		qf = self.ismall2qsmall(ifn)
		qr = self.ismall2qsmall(irn)

		Ids = CS_FACTOR* self.NPMOS * device.L/Leff *self.scaling.Is * (ifn - irn)
		
		vd_real = vd if CS_FACTOR == 1 else vs
		vs_real = vs if CS_FACTOR == 1 else vd

		opdict.update({'state':(vd_real*self.NPMOS, vg*self.NPMOS, vs_real*self.NPMOS)})
		opdict.update({'Ids':Ids, "Weff":device.W, "Leff":Leff, 'Vp':VP})
		opdict.update({'ifn':ifn, "irn":irn, "n":nv, 'beta':.5*self.KP*device.W/Leff, 'Ispec':self.scaling.Is})
		opdict.update({'VTH':self.VTO, "VOD":self.NPMOS*nv*(VP-VS), 'SAT':ifn>irn*self.SATLIM})
		opdict.update({'qf':qf*self.scaling.Qs, 'qr':qr*self.scaling.Qs})

		if Ids > self.WMSI_factor*self.scaling.Is: 
			WMSI = 2
		elif Ids < self.scaling.Is/self.WMSI_factor:
			WMSI = 0
		else:
			WMSI = 1
		opdict.update({'WMSI':WMSI})

		if debug: print "current:", Ids
		
		return Ids, qf, qr

	def get_gms(self, device, (vd, vg, vs), opdict=None, debug=False):
		(j1, j2, j3), CS_FACTOR = self.get_voltages(vd, vg, vs)
		Ids, qf, qr = self.get_ids(device, (vd, vg, vs), opdict, debug)
		if CS_FACTOR == +1:
			gms = -1.0*self.scaling.Gs*qf
		elif CS_FACTOR == -1:
			gms = self.scaling.Gs*qr
		return gms

	def get_gmd(self, device, (vd, vg, vs), opdict=None, debug=False):
		(j1, j2, j3), CS_FACTOR = self.get_voltages(vd, vg, vs)
		Ids, qf, qr = self.get_ids(device, (vd, vg, vs), opdict, debug)
		if CS_FACTOR == +1:
			gmd = self.scaling.Gs*qr
		elif CS_FACTOR == -1:
			gmd = self.scaling.Gs*qf
		return gmd

	def get_gmg(self, device, (vd, vg, vs), opdict=None, debug=False):
		VP, nv, nq = self.get_vp_nv_nq(float(vg))
		Ids, qf, qr = self.get_ids(device, (vd, vg, vs), opdict, debug)
		(j1, j2, j3), CS_FACTOR = self.get_voltages(vd, vg, vs)
		gmg = CS_FACTOR*self.scaling.Gs*(qf-qr)/nv
		return gmg

	def get_ismall(self, vsmall, ip_abs_err, iguess=None, debug=False):
		if iguess is None:
			iguess = 1
		if math.isnan(vsmall):
			raise Exception, "vsmall is NaN!!"
		check = False
		ismall = iguess
		if debug: iter_c = 0
		if debug: print "get_ismall IN"
		#if vsmall < utilities.EPS:
		#	print "vsmall:", vsmall, "returning EPS"
		#	return utilities.EPS
			
		while True:
			if debug: iter_c = iter_c + 1
			vsmall_iter = self.get_vsmall(ismall)
			if debug: print " ->", ismall, vsmall
			deltai = (vsmall - vsmall_iter)/self.get_dvsmall_dismall(ismall)
			#if ismall == 0:
			#	ismall = utilities.EPS
			#	if debug: print "ismall = 0!"
			if (abs(deltai) < ip_abs_err or abs(deltai) <= abs(ismall)*options.ier):# or \
				if not check:
					check = True
				else:
					break
			else:
				check = False
			if deltai == 0:
				break
			if math.isnan(ismall):
				print "Ismall is NaN!!"
				exit()
			if ismall == 0:
				check = True
			else:
				ratio = deltai/ismall	
				if ratio > 3: 			 
					ismall = 2*ismall
				elif ratio <= -1:
					ismall = 0.1*ismall
				else:			
					ismall = ismall + deltai
		if debug: 
			print str(iter_c) + " iterations."
			print ismall, ismall/iguess
		if debug: print "get_ismall OUT"
		assert ismall >= 0
		return ismall
	
	def get_vsmall(self, ismall, verbose=3):
		#print "AA ", 0.25 + ismall, -0.5 + math.sqrt(0.25 + ismall)
		if abs(ismall) < utilities.EPS:
			ismall = (1 - 2*(ismall<0))*utilities.EPS
			if verbose == 6:
				print "EKV: Machine precision limited the resolution on i. (i<EPS)"

		vsmall = math.log(math.sqrt(.25 + ismall) - 0.5) + 2*math.sqrt(.25 + ismall) - 1.0
		return vsmall

	def get_dvsmall_dismall(self, ismall, verbose=3):
		if abs(ismall) < utilities.EPS:
			ismall = (1-2*(ismall<0))*utilities.EPS
			if verbose == 6:
				print "EKV: Machine precision limited the resolution on dv/di in the NR iteration. (i<EPS)"
		dvdi = 1.0/(math.sqrt(.25+ismall)-.5) * .5/math.sqrt(.25 + ismall) + 1.0/math.sqrt(.25 + ismall)
		return dvdi

	def ismall2qsmall(self, ismall, verbose=0):
		if verbose == 6:
			print "EKV: Machine precision limited the resolution on q(s,d). (i<EPS)"
		qsmall = math.sqrt(.25+ismall) - .5
		return qsmall

	def qsmall2ismall(self, qsmall):
		ismall = qsmall**2 + qsmall 
		return ismall

	def _self_check(self):
		ret = True, ""

		if self.NSUB is not None and self.NSUB < 0:
			ret = (False, "NSUB")
		elif self.U0 is not None and not self.U0 > 0:
			ret = (False, "UO")
		elif not self.GAMMA > 0:			 
			ret = (False, "GAMMA")
		elif not self.PHI > 0.1:			 
			ret = (False, "PHI")
		return ret

	def _device_check(self, adev):
		if not adev.L > 0:
			ret = (False, "L")
		elif not adev.W > 0:
			ret = (False, "W")
		# we don't check for fractional series/shunt devices
		# DIY 
		elif not adev.N > 0:
			ret = (False, "N")
		elif not adev.M > 0:
			ret = (False, "M")
		else:
			ret = (True, "")
		return ret
		
if __name__ == '__main__':
	import matplotlib.pyplot as plt
	ekv_m = ekv_mos_model(TYPE='n', KP=50e-6, VTO=.4)
	ma = ekv_device(1, 2, 3, 4, W=10e-6,L=1e-6, model=ekv_m)
	ma.descr = "1"

	vd = 0
	vg = 1
	vs = 0
	ma.print_op_info(((vd, vg, vs),))

	if True:
		data0 = []
		data1 = []
		vs = 2.5
		if True:
			vd = 1
			for Vhel in range(250):
				vg = (Vhel+1)/100.0
				ma.update_status_dictionary(((vd, vg, 0),))						
				data0.append(ma.opdict['Ids'])
				#print "Current for vd", vd, "vg", vg, "vs", vs
				data1.append(ma.opdict['TEF'])
		plt.plot(data0, data1)
		plt.title('EKV MODEL CHECK: VD_SWEEP')
		plt.legend(['GM/ID'])
		plt.show()



	#print data2
