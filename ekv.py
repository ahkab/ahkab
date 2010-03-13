# -*- coding: iso-8859-1 -*-
# ekv.py
# EKV MOS transistor model
# Copyright 2010 Giuseppe Venturini
# 
# The EKV model was developed by Matthias Bucher, Christophe Lallement, 
# Christian Enz, Fabien Théodoloz, François Krummenacher at the Electronics 
# Laboratories, Swiss Federal Institute of Technology (EPFL), 
# Lausanne, Switzerland
# and the Tecnical Report upon which this implementation is based is here:
# <http://legwww.epfl.ch/ekv/pdf/ekv_v262.pdf>.
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
Lausanne, Switzerland
and the Tecnical Report upon which this implementation is based is here:
<http://legwww.epfl.ch/ekv/pdf/ekv_v262.pdf>.

This module defines two classes:
 ekv_device
 ekv_mos_model

and several default values for the model. Their meaning is stated in the 
model doc.

TODO:
- Complex mobility degradation is missing
- Transcapacitances
- Quasistatic implementation
"""

import constants, options, utilities, printing
import math, random

# DEFAULT VALUES FOR 500n CH LENGTH
COX_DEFAULT = .7e-3
XJ_DEFAULT = .1e-6

DW_DEFAULT = 0
DL_DEFAULT = 0
VTO_DEFAULT = .5
GAMMA_DEFAULT = 1
PHI_DEFAULT = .7
KP_DEFAULT = 50e-6
E0_DEFAULT = 1e12
UCRIT_DEFAULT = 2e6

THETA_DEFAULT = 0
LAMBDA_DEFAULT = .5
WETA_DEFAULT = 0#.25
LETA_DEFAULT = 0#.1
Q0_DEFAULT = 0
LK_DEFAULT = .29e-6

IBA_DEFAULT = 0
IBB_DEFAULT = 3e8
IBN_DEFAULT = 1

TCV_DEFAULT = 1e-3
BEX_DEFAULT = -1.5 
UCEX_DEFAULT = .8
IBBT_DEFAULT = 9e-4



class ekv_device:
	# enable this to use the pseudo-random generator
	OPT_SPICE_IT_UP = False
	INIT_IFRN_GUESS = 1e-3
	def __init__(self, nd, ng, ns, nb, L, W, model, M=1, N=1):
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
		self.device.M = M #parallel multiple device number
		self.device.N = N #series multiple device number
		self.ekv_model = model
		self.opdict = {}
		self.opdict.update({'state':(float('nan'), float('nan'), float('nan'))})
		self.opdict.update({'ifn':self.INIT_IFRN_GUESS})
		self.opdict.update({'irn':self.INIT_IFRN_GUESS})
		self.opdict.update({'ip_abs_err':self.ekv_model.get_ip_abs_err(self.device)})
		self.letter_id = 'M'
		self.is_nonlinear = True
		self.dc_guess = [self.ekv_model.VTO*(0.1)*self.ekv_model.NPMOS, self.ekv_model.VTO*(1.1)*self.ekv_model.NPMOS, 0]
		self.randomstate = None
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
		return ((self.n1, self.n2),(self.n1, self.nb))
	
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
		if op_index == 0: 		
			ret, j1, j2 = \
				self.ekv_model.get_ids(self.device, ports_v, \
				self.opdict)
		elif op_index == 1:
			ret, j1, j2 = \
				self.ekv_model.get_idb(self.device, ports_v, \
				self.opdict)
		return ret

	def print_op_info(self, ports_v):
		"""Operating point info, for design/verification. """
		mos_type = self._get_mos_type()
		
		gmd, gmg, gms, ids = (None, None, None, None)
		if not (self.opdict['state'] == ports_v[0] and self.opdict.has_key('gmd')):
			gmd = self.g(1, ports_v[0], 0) + self.g(0, ports_v[0], 0)
		else:	
			gmd = self.opdict['gmd'] + self.g(1, ports_v, 0)		
		if not (self.opdict['state'] == ports_v[0] and self.opdict.has_key('gmg')):
			gmg = self.g(0, ports_v[0], 1)
		else:	
			gmg = self.opdict['gmg']		
		if not (self.opdict['state'] == ports_v[0] and self.opdict.has_key('gms')):
			gms = self.g(0, ports_v[0], 2)
		else:	
			gms = self.opdict['gms']		
		if not (self.opdict['state'] == ports_v[0] and self.opdict.has_key('Ids')):
			ids = self.i(0, ports_v[0])
		else:	
			ids = self.opdict['Ids']
		if not (self.opdict['state'] == ports_v[0] and self.opdict.has_key('Idb')):
			idb = self.i(1, ports_v[0])
		else:	
			idb = self.opdict['Idb']

		if ids == 0:
			TEF = float('nan')
		else:
			TEF = gms*constants.Vth()/ids	
		
		status = "SATURATION" if self.opdict['SAT'] else "LINEAR MODE"

		arr = [["M"+self.descr, mos_type.upper()+" ch",status, "", "", "", "", "", "", "", "",""],]
		arr.append(["Weff", "[m]:", str(self.opdict['Weff'])+" ("+str(self.device.W)+")", "Leff", "[m]:", str(self.opdict['Leff'])+ " ("+str(self.device.L)+")", "M:", "", self.device.M, "N:", "", self.device.N])
		arr.append(["Vdb", "[V]:", float(ports_v[0][0]), "Vgb", "[V]:", float(ports_v[0][1]), "Vsb", "[V]:", float(ports_v[0][2]),  "Vp", "[V]:", self.opdict['Vp'],])
		arr.append([ "VTH", "[V]:", self.opdict['VTH'], "VOD", "[V]:", self.opdict['VOD'], "VDSAT", "[V]:", self.opdict['VDSAT'], "VA", "[V]:", str(ids/gmd)])
		arr.append(["Ids", "[A]:", ids, "Idb", "[A]:", idb, "", "", "", "TEF:", "", str(TEF),]) 
		arr.append(["gmg", "[S]:", gmg, "gms", "[S]:", gms, "rob", "[Ohm]:", 1/gmd, "", "", ""])
		arr.append(["if:", "", self.opdict['if'],"ir:", "", self.opdict['ir'], "n: ", "",self.opdict['n'],"beta:", "", self.opdict['beta']])
		#arr.append([  "", "", "", "", "", ""])

		printing.table_print(arr)

#		print "M"+self.descr+":", mos_type.upper(), status
#		print "  ", "Weff:", self.opdict['Weff'], "("+str(self.device.W)+")", "Leff:", self.opdict['Leff'], "("+str(self.device.L)+")"
#		print "  ", "Vdb:", str(ports_v[0][0]), "vgb:", str(ports_v[0][1]), "vsb:", str(ports_v[0][2])
#		print "  ", 'Vp', self.opdict['Vp'], 'VTH', self.opdict['VTH'], 'VOD', self.opdict['VOD'], 'VDSAT', self.opdict['VDSAT']
#		print "  ", "ids =", str(ids), "idb =", str(idb) 
#		print "  ", "gmg:", str(gmg), "gms:", str(gms), "rob:", str(1/gmd)
#		print "  ", "(if:", str(self.opdict['if'])+ ", ir:", str(self.opdict['if'])+", n: "+str(self.opdict['n'])+", beta: "+str(self.opdict['beta'])+")"
#		print "  ", "TEF =", str(TEF), "VA:", str(ids/gmd)		

	
	def g(self, op_index, ports_v, port_index, time=0):
		"""Returns the differential (trans)conductance rs the port specified by port_index
		when the element has the voltages specified in ports_v across its ports,
		at (simulation) time.

		This is here not computed through analytical derivation rather it is
		approximated by [I2(v+dv)-I1(v)]/dv
		
		ports_v: a list in the form: [voltage_across_port0, voltage_across_port1, ...]
		port_index: an integer, 0 <= port_index < len(self.get_ports())
		time: the simulation time at which the evaluation is performed. Set it to
		None during DC analysis.
		"""
		if self.randomstate is not None and self.OPT_SPICE_IT_UP:			
			random.setstate(self.randomstate)
		step = 	math.sqrt(options.ver)*(1 + self.OPT_SPICE_IT_UP*random.random())#*base_step
		if self.OPT_SPICE_IT_UP: 
			self.randomstate = random.getstate()		
				
		vd_2 = ports_v[0] if port_index != 0 else ports_v[0]+step
		vg_2 = ports_v[1] if port_index != 1 else ports_v[1]+step
		vs_2 = ports_v[2] if port_index != 2 else ports_v[2]+step
		vd_1 = ports_v[0]
		vg_1 = ports_v[1]
		vs_1 = ports_v[2]

		ports_v2 = (vd_2, vg_2, vs_2)
		ports_v1 = (vd_1, vg_1, vs_1)

		if op_index == 0:		
			I2, j1, j2 = self.ekv_model.get_ids(self.device, ports_v2, self.opdict)
			I1, j1, j2 = self.ekv_model.get_ids(self.device, ports_v1, self.opdict)
		else:
			I2, j1, j2 = self.ekv_model.get_idb(self.device, ports_v2, self.opdict)
			I1, j1, j2 = self.ekv_model.get_idb(self.device, ports_v1, self.opdict)

		g = (I2-I1)/step

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

class ekv_mos_model:
	Cepsi = 4*(22e-3)**2
	Ca = 0.028	
	def __init__(self, name=None, TYPE='n', TNOM=None, COX=None, XJ=None, DW=None, DL=None,\
	GAMMA=None, NSUB=None, PHI=None, VTO=None, KP=None, E0=None, UCRIT=None, \
	TOX=None, VFB=None, U0=None, VMAX=None, THETA=None, LAMBDA=None, \
	WETA=None, LETA=None, Q0=None, LK=None, IBA=None, IBB=None, IBN=None, \
	TCV=None, BEX=None, UCEX=None, IBBT=None):
		
		self.name = "model_ekv0" if name is None else name
		Vth = constants.Vth()
		self.TNOM = TNOM if TNOM is not None else constants.Tref	
		#print "TYPE IS:" + TYPE		
		self.NPMOS = 1 if TYPE == 'n' else -1

		# optional parameters (no defaults)	
		self.TOX = TOX
		self.NSUB = NSUB
		self.VFB = self.NPMOS*VFB if VFB is not None else None
		self.U0 = U0
		self.VMAX = VMAX
		self.THETA = THETA if THETA is not None else THETA_DEFAULT

		# crucial parameters
		if COX is not None: 
			self.COX = COX
		elif TOX is not None:
			self.COX = constants.si.esi/TOX
		else:
			self.COX = COX_DEFAULT
		self.XJ = XJ if XJ is not None else XJ_DEFAULT

		self.DW = DW if DW is not None else DW_DEFAULT
		self.DL = DL if DL is not None else DL_DEFAULT
	
		if GAMMA is not None:
			self.GAMMA = GAMMA
		elif NSUB is not None:
			math.sqrt(2*constants.e*constants.si.esi*NSUB*10**6/self.COX)
		else:
			self.GAMMA = GAMMA_DEFAULT
		if PHI is not None:
			self.PHI = PHI
		elif NSUB is not None:
			self.PHI = 2*constants.Vth(self.TNOM)*math.log(NSUB*10**6/constants.si.ni(self.TNOM))
		else:
			self.PHI = PHI_DEFAULT
		if VTO is not None:
			self.VTO = self.NPMOS*VTO
		elif VFB is not None:
			self.VTO = VFB + PHI + GAMMA*PHI #inv here??
		else:
			self.VTO = self.NPMOS*VTO_DEFAULT

		if KP is not None:
			self.KP = KP
		elif U0 is not None:
			self.KP = (U0*10**-4)*self.COX
		else:
			self.KP = KP_DEFAULT
	
		if E0 is not None:
			self.E0 = E0 
		elif THETA is not None: 
			self.E0 = 0
		else:
			self.E0 = E0_DEFAULT
	
		if UCRIT is not None:
			self.UCRIT = UCRIT 
		elif VMAX is not None and U0 is not None:
			self.UCRIT
		else:
			self.UCRIT = UCRIT_DEFAULT 
	
		#Channel length modulation and charge sharing parameters
		self.LAMBDA = LAMBDA if LAMBDA is not None else LAMBDA_DEFAULT
		self.WETA = WETA if WETA is not None else WETA_DEFAULT
		self.LETA = LETA if LETA is not None else LETA_DEFAULT
		# Reverse short-channel effect parameters
		self.Q0 = Q0 if Q0 is not None else Q0_DEFAULT
		self.LK = LK if LK is not None else LK_DEFAULT
		# impact ionization current parameters
		self.IBA = IBA if IBA is not None else IBA_DEFAULT
		self.IBB = IBB if IBB is not None else IBB_DEFAULT
		self.IBN = IBN if IBN is not None else IBN_DEFAULT
		# Intrinsic model temperature parameters
		self.TCV = self.NPMOS*TCV if TCV is not None else self.NPMOS*TCV_DEFAULT
		self.BEX = BEX if BEX is not None else BEX_DEFAULT
		self.UCEX = UCEX if UCEX is not None else UCEX_DEFAULT
		self.IBBT = IBBT if IBBT is not None else IBBT_DEFAULT
	
		# temperature effects update
		self.VTO = self.VTO - self.TCV*(constants.T-self.TNOM)
		self.KP = self.KP*(constants.T/self.TNOM)**self.BEX
		self.UCRIT = self.UCRIT*(constants.T/self.TNOM)**self.UCEX
		self.PHI = self.PHI * constants.T/self.TNOM + 3*constants.Vth(self.TNOM)*math.log(constants.T/self.TNOM) \
			   - constants.si.Eg(self.TNOM)*constants.T/self.TNOM + constants.si.Eg(constants.T)
		self.IBB = self.IBB*(1+ self.IBBT*(constants.T - self.TNOM))
	
		#Setup switches
		self.XQC = 1
		self.SATLIM = math.exp(4)
		self.NQS = 0

		sc, sc_reason = self._self_check()
		if not sc:
			raise Exception, sc_reason + " out of range"

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
		return ((float(vd_new), float(vg), float(vs_new)), cs)

	def get_ip_abs_err(self, device):
		return  options.iea / (4*constants.Vth()**2*self.KP*device.M*device.W/device.L)

	def _get_Vgp(self, vg, Weff, Leff, debug=False):
		"""Do not forget to setup the inputs FIRST with get_voltages!
		"""
		# Reverse short-channel effect (RSCE)
		Psi = self.Ca * (10*Leff/self.LK - 1) # eq. 31
		Delta_vrsce = 2*self.Q0/self.COX/(1 + .5*(Psi + math.sqrt(Psi**2+self.Cepsi)))**2 # eq. 32
		if debug: 
			print "Psi:", Psi, "Delta_vrsce:", Delta_vrsce

		# Effective gate voltage including RSCE
		Vgp = vg -self.VTO -Delta_vrsce +self.PHI +self.GAMMA*math.sqrt(Psi) # eq. 33

		return Vgp, Delta_vrsce
	
	def get_ids(self, device, (vd, vg, vs), opdict=None, debug=False):
		# internal variables setup	
		if debug: print "Current for vd:", vd, "vg:", vg, "vs:", vs
		ip_abs_err = self.get_ip_abs_err() if opdict['ip_abs_err'] is None else opdict['ip_abs_err']
		Vth = constants.Vth()
		(vd, vg, vs), CS_FACTOR = self.get_voltages(vd, vg, vs)	
		vds = .5*(vd - vs) # eq. 49		
		Weff, Leff = self.get_eff_wl(device.W, device.L)
		
		if debug:
			print "vd:", vd, "vg:",vg, "vs:", vs, "vds:", vds
			print "ip_abs_err:", ip_abs_err
			print "Vth:", Vth
			print "Weff:", Weff, "Leff:", Leff

		Vgp, Delta_vrsce = self._get_Vgp(vg, Weff, Leff)

		# Effective substrate factor including charge-sharing for short and narrow channels
		if Vgp > 0: 
			Vp0 = Vgp -self.PHI -self.GAMMA*(math.sqrt(Vgp+(self.GAMMA/2)**2)-self.GAMMA/2)
		else:
			Vp0 = -self.PHI
		Vsp = .5 * (vs + self.PHI + math.sqrt((vs + self.PHI)**2+(4*Vth)**2))
		Vdp = .5 * (vd + self.PHI + math.sqrt((vd + self.PHI)**2+(4*Vth)**2))
		
		if debug:
			print "Vsp:",Vsp, "Vdp:", Vdp, "CK:", math.sqrt(self.PHI + Vp0)
		
		gam0 = self.GAMMA
		if self.LETA > 0:
			gam0 = gam0 - constants.si.esi/self.COX*self.LETA/Leff*(math.sqrt(Vsp)+math.sqrt(Vdp))
			if debug: print gam0
		if self.WETA > 0:
			gam0 = gam0 - constants.si.esi/self.COX*3*self.WETA/Weff*math.sqrt(self.PHI + Vp0)
		gamp = .5 * (gam0 + math.sqrt(gam0**2+0.1*Vth)) #FIXME
		if debug: print "gam0:", gam0, "gamp:", gamp

		# Pinch-off voltage including short and narrow channel effects, eq. 38
		if Vgp > 0:
			Vp = Vgp - self.PHI - gamp*(math.sqrt(Vgp + (gamp/2)**2)-gamp/2)
		else:
			Vp = -self.PHI

		assert not math.isnan(Vp) 

		# slope factor, eq. 39		
		n = 1 + self.GAMMA/(2*math.sqrt(Vp +self.PHI+4*Vth))

		# forward current scaled control voltage
		v_if = (Vp - vs)/Vth
		#assert v_if > 0
		# scaled forward current
		ifn = self.get_ismall(v_if, opdict['ip_abs_err'], opdict['ifn'])
		if debug: print "# ifn OK"		

		if ifn != 0:		
			# velocity saturation voltage, eq. 45
			Vc = self.UCRIT*device.N*Leff

			vdss  = Vc * (math.sqrt(.25+Vth/Vc *math.sqrt(ifn))-.5) # eq. 46
			# Drain-to-source saturation voltage for reverse normalized current, eq. 47		
			vdssp = Vc * (math.sqrt(.25 +Vth/Vc *(math.sqrt(ifn)-.75*math.log(ifn)))-.5) + Vth*(math.log(Vc/(2*Vth))-.6)

			# channel length modulation
			assert not math.isnan(vdss)
			#print math.sqrt(ifn), vdss/Vth
			#assert math.sqrt(ifn) > vdss/Vth
			vser_1 = math.sqrt(ifn) - vdss/Vth
			if vser_1 < 0: vser_1 = 0
			delta_v = 4*Vth*math.sqrt(self.LAMBDA*(vser_1)+1/64) # eq. 48
			Vip = math.sqrt(vdss**2 + delta_v**2) - math.sqrt((vds - vdss)**2 + delta_v**2) #eq 50
			Lc = math.sqrt(constants.si.esi*self.XJ/self.COX) #eq. 51
			delta_l = self.LAMBDA * Lc * math.log(1+ (vds - Vip)/(Lc*self.UCRIT)) #eq. 52
	
			# Equivalent channel length including channel-length modulation and velocity saturation
			Lp = device.N*Leff - delta_l + (vds + Vip)/self.UCRIT #eq. 53
			Lmin = device.N*Leff/10 #eq. 54
			Leq = .5*(Lp + math.sqrt(Lp**2 + Lmin**2)) #eq. 55
			assert not math.isnan(vdssp)
			assert not math.isnan(delta_v)
			if debug: print "# channel length modulation OK"

			v_irp = (Vp - vds + vs - math.sqrt(vdssp**2 + delta_v**2)+math.sqrt((vds-vdssp)**2+delta_v**2))/Vth
			if debug: print delta_v

			irn = self.get_ismall(v_irp, opdict['ip_abs_err'], opdict['irn'])
		else:
			vdss = 0
			irn = opdict['ip_abs_err']*1.1
			ifn = opdict['ip_abs_err']*1.1*2
			Leq = Leff
			
		beta0 = self.KP*device.M*Weff/Leq
		
		if self.THETA == 1:
			Vpp = .5 * (Vp + math.sqrt(Vp**2 + 2.0*Vth**2))
			beta = beta0/(1+self.THETA*Vpp)
		else:
			beta = beta0
			#FIXME Complex mobility degradation missing

		Is = 2* n * beta * Vth**2
		Ids = CS_FACTOR* self.NPMOS * Is*(ifn - irn)
		
		vd_real = vd if CS_FACTOR == 1 else vs
		vs_real = vs if CS_FACTOR == 1 else vd
		opdict.update({'state':(vd_real*self.NPMOS, vg*self.NPMOS, vs_real*self.NPMOS)})
		opdict.update({'Ids':Ids, "Weff":Weff, "Leff":Leff, 'Vp':Vp})
		opdict.update({'if':ifn, "ir":irn, "n":n, 'beta':beta})
		opdict.update({'VTH':self.VTO + Delta_vrsce + gamp*math.sqrt(Vsp)-self.GAMMA*math.sqrt(self.PHI), \
				"VOD":self.NPMOS*n*(Vp-vs), "VDSAT":2*vdss, 'SAT':ifn>irn*self.SATLIM})
		if debug: print "current:", Ids
		return Ids, vdss, Lc

	def get_idb(self, device, (vd, vg, vs), opdict, debug=False):
		ids, vdss, Lc = self.get_ids(device, (vd, vg, vs), opdict, debug)
		Vib = vd - vs - self.IBN * 2 * vdss
		if Vib > 0 and self.IBA != 0: 
			idb = ids * self.IBB/self.IBA * Vib * math.exp(-self.IBB*Lc/Vib)
		else:
			idb = 0	
		opdict.update({'Idb':idb})
		return idb, vdss, Lc
		

	def get_eff_wl(self, w, l):
		return w+self.DW, l+self.DL


	def get_ismall(self, vsmall, ip_abs_err, iguess=None, debug=False):
		if iguess is None:
			iguess = 1
		if math.isnan(vsmall):
			raise Exception, "vsmall is NaN!!"
		check = False
		ismall = iguess
		#if debug: iter_c = 0

		#if vsmall < utilities.EPS:
		#	print "vsmall:", vsmall, "returning EPS"
		#	return utilities.EPS
			
		while True:
			#if debug: iter_c = iter_c + 1			
			vsmall_iter = self.get_vsmall(ismall)
			if debug: print vsmall_iter, vsmall
			deltai = (vsmall - vsmall_iter)/self.get_dvsmall_dismall(ismall)
			ratio = deltai/ismall			
			if ismall == 0:
				ismall = utilities.EPS
			#print ismall, deltai
			elif abs(deltai) < ip_abs_err or abs(deltai/ismall) < options.ier:
				if not check:
					check = True
				else:
					break
			else:
				check = False
			if math.isnan(ismall):
				print "Ismall is NaN!!"
				exit()
			if ratio > 3: 			 
				ismall = 3*deltai
			elif ratio <= -1:
				ismall = 0.1*ismall
			else:			
				ismall = ismall + deltai
		#if debug: print str(iter_c) + " iterations."
		return ismall
	
	def get_vsmall(self, ismall, verbose=3):
		#print "AA ", 0.25 + ismall, -0.5 + math.sqrt(0.25 + ismall)
		if abs(ismall) < utilities.EPS:
			ismall = (1-2*(ismall<0))*utilities.EPS
			if verbose == 6:
				print "EKV: Machine precision problem detected in current evaluation"
		vsmall = -1.0 + math.log(-0.5 + .5*math.sqrt(1 + 4*ismall)) + math.sqrt(1 + 4*ismall)
		return vsmall

	def get_dvsmall_dismall(self, ismall, verbose=3):
		if abs(ismall) < utilities.EPS:
			ismall = (1-2*(ismall<0))*utilities.EPS
			if verbose == 6:
				print "EKV: Machine precision problem detected in the NR iteration"
		dvdi = -1.0/(.5*math.sqrt(1 + 4*ismall)*(1.0 - math.sqrt(1 + 4*ismall))) + 2/math.sqrt(1 + 4*ismall)
		return dvdi

	def _self_check(self):
		if not self.XJ > 0:
			ret = (False, "XJ")
		elif self.NSUB is not None and not self.NSUB > 0:
			ret = (False, "NSUB")
		elif self.U0 is not None and not self.U0 > 0:
			ret = (False, "UO")
		elif self.VMAX is not None and not self.VMAX > 0:
			ret = (False, "VMAX")
		elif not self.THETA >= 0:			 
			ret = (False, "THETA")
		elif not self.GAMMA > 0:			 
			ret = (False, "GAMMA")
		elif not self.PHI > 0.1:			 
			ret = (False, "PHI")
		elif not self.E0 > 1e5:			 
			ret = (False, "E0")
		elif not self.UCRIT > 0:			 
			ret = (False, "UCRIT")
		elif not self.IBN > 0.1:			 
			ret = (False, "IBN")
		elif not self.IBB > 1e8:			 
			ret = (False, "IBB")
		elif not self.LK > 0:			 
			ret = (False, "LK")
		else:
			ret = True, ""
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
	ekv_m = ekv_mos_model(TYPE='n')
	ma = ekv_device(1, 2, 3, 4, 10e-6,1e-6, ekv_m)
	ma.descr = "1"

	print "#Vg\tId\tgmd\tgmg\tgms"
	vdb = 1
	for Vhel in range(250):
		Vg = Vhel/100.0
		str_mine =  str(Vg)+"\t"+str(ma.i(0, (vdb, Vg, 0)))+"\t"+str(ma.g(0, (vdb, Vg, 0), 0))+"\t"+str(ma.g(0, (vdb, Vg, 0), 1))+"\t"+str(ma.g(0, (vdb, Vg, 0), 2))
		print str_mine
	print ma.print_op_info(((0.33031697099518587, 1.1704486464618264, 0.0),(0.33031697099518587, 1.1704486464618264, 0.0)))
