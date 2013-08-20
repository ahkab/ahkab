# -*- coding: iso-8859-1 -*-
# diode.py
# Diode model
# Copyright 2006-2013 Giuseppe Venturini

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

import math, numpy
import constants, printing, dc_analysis, utilities

"""
Contains a diode element and its model class.
                                  
         |\|                        
  n1 o---| |---o n2
         |/|                        
                                  
                                  
"""

class diode:
	letter_id = "d"
	is_nonlinear = True
	is_symbolic = True
	dc_guess = [0.425]
	def __init__(self, n1, n2, model, AREA=None, T=None, ic=None, off=False):
		class dev_class: pass
		self.device = dev_class()
		self.device.AREA = AREA if AREA is not None else 1.0
		self.device.T = T
		self.device.last_vd = .425
		self.n1 = n1
		self.n2 = n2
		self.ports = ((self.n1, self.n2),)
		self.model = model
		if self.device.T is None:
			self.device.T = constants.T
		
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
		return self.device.T
	def set_temperature(self, T):
		"""Set the operating temperature IN KELVIN degrees"
		self.device.T = T
	def __str__(self):
		T = self._get_T()
		rep = "%s area=%g T=%g" % (self.model.name, self.device.AREA, self.device.T)
		if self.ic is not None:
			rep = rep + " ic="+str(self.ic)
		elif self.off:
			rep += " off"
		return rep

	def get_output_ports(self):
		return self.ports

	def get_drive_ports(self, op):
		return self.ports

	def i(self, op_index, ports_v, time=0): #with gmin added
		v = ports_v[0]
		i = self.model.get_i(op_index, ports_v, self.device)
		return i

	def g(self, op_index, ports_v, port_index, time=0):
		if not port_index == 0: 
			raise Exception, "Attepted to evaluate a diode's gm on an unknown port."
		return self.model.get_gm(op_index, ports_v, port_index, self.device)

	def get_op_info(self, ports_v_v):
		vn1n2 = float(ports_v_v[0][0])
		idiode = self.i(0, (vn1n2,))
		gmdiode = self.g(0, (vn1n2,), 0)
		info = ["V(n1-n2): ", vn1n2, "[V]", "I(n1-n2):", idiode, "[A]", "P:", vn1n2*idiode, "g:", gmdiode, "[A/V]", "T:", self._get_T(), "K" ]
		arr = [[self.letter_id.upper()+self.descr] + info]
		strarr = printing.table_setup(arr)
		return strarr

	def print_op_info(self, ports_v):
		print self.get_op_info(ports_v),

IS_DEFAULT = 1e-14 # A
N_DEFAULT = 1.0
ISR_DEFAULT = 0.0 #A
NR_DEFAULT = 2.0
RS_DEFAULT = 0.0 # ohm
CJ0_DEFAULT = 0.0 
M_DEFAULT = .5
VJ_DEFAULT = .7
FC_DEFAULT = .5
CP_DEFAULT = .0
TT_DEFAULT = .0
BV_DEFAULT = float('inf')
IBV_DEFAULT = 1e-3
KF_DEFAULT = .0
AF_DEFAULT = 1.
FFE_DEFAULT = 1.
TEMP_DEFAULT = 26.85
XTI_DEFAULT = 3.0
EG_DEFAULT = 1.11
TBV_DEFAULT = 0.0
TRS_DEFAULT = 0.0
TTT1_DEFAULT = 0.0
TTT2_DEFAULT = 0.0
TM1_DEFAULT = 0.0
TM2_DEFAULT = 0.0
T_DEFAULT = utilities.Celsius2Kelvin(26.85)
AREA_DEFAULT = 1.0
class diode_model:
	def __init__(self, name, IS=None, N=None, ISR=None, NR=None, RS=None, 
		     CJ0=None, M=None, VJ=None, FC=None, CP=None, TT=None, 
		     BV=None, IBV=None, KF=None, AF=None, FFE=None, TEMP=None, 
		     XTI=None, EG=None, TBV=None, TRS=None, TTT1=None, TTT2=None, 
		     TM1=None, TM2=None):
		self.name = name
		self.IS = float(IS) if IS is not None else IS_DEFAULT
		self.N = float(N) if N is not None else N_DEFAULT
		self.ISR = float(ISR) if ISR is not None else ISR_DEFAULT
		self.NR = float(NR) if NR is not None else NR_DEFAULT
		self.RS = float(RS) if RS is not None else RS_DEFAULT
		self.CJ0 = float(CJ0) if CJ0 is not None else CJ0_DEFAULT
		self.M = float(M) if M is not None else M_DEFAULT
		self.VJ = float(VJ) if VJ is not None else VJ_DEFAULT
		self.FC = float(FC) if FC is not None else FC_DEFAULT
		self.CP = float(CP) if CP is not None else CP_DEFAULT
		self.TT = float(TT) if TT is not None else TT_DEFAULT
		self.BV = float(BV) if BV is not None else BV_DEFAULT
		self.IBV = float(IBV) if IBV is not None else IBV_DEFAULT
		self.KF = float(KF) if KF is not None else KF_DEFAULT
		self.AF = float(AF) if AF is not None else AF_DEFAULT
		self.FFE = float(FFE) if FFE is not None else FFE_DEFAULT
		self.TEMP = utilities.Celsius2Kelvin(float(TEMP)) if TEMP is not None else TEMP_DEFAULT
		self.XTI = float(XTI) if XTI is not None else XTI_DEFAULT
		self.EG = float(EG) if EG is not None else EG_DEFAULT
		self.TBV = float(TBV) if TBV is not None else TBV_DEFAULT
		self.TRS = float(TRS) if TRS is not None else TRS_DEFAULT
		self.TTT1 = float(TTT1) if TTT1 is not None else TTT1_DEFAULT
		self.TTT2 = float(TTT2) if TTT2 is not None else TTT2_DEFAULT
		self.TM1 = float(TM1) if TM1 is not None else TM1_DEFAULT
		self.TM2 = float(TM2) if TM2 is not None else TM2_DEFAULT
		self.T = T_DEFAULT
		self.last_vd = None
		self.VT = constants.Vth(self.T)
	def get_i(self, op_index, ports_v, dev):
		if dev.T != self.T:
			self.set_temperature(dev.T)
		if not self.RS:
			i = self._get_i(ports_v[0])*dev.AREA if self.RS == 0 else self._get_irs(ports_v, dev)
			dev.last_vd = ports_v[0]
		else:
			i = self._get_irs(ports_v, dev)
		return i
	def _get_irs(self, ports_v, dev):
		vth = self.VT
		vd = dev.last_vd if dev.last_vd is not None else vth
		RS = self.RS/dev.AREA
		idiode = self._get_i(vd)*dev.AREA
		while True:	
			gm = self.get_gm(0, [vd], 0, dev, rs=False)
			dvd = (ports_v[0] - idiode*self.RS-vd)/(1.0 + gm*RS)
			vd = vd + min(self.VT, abs(dvd))*numpy.sign([dvd])[0]
			idiode_old = idiode
			idiode = self._get_i(vd)*dev.AREA
			di = idiode - idiode_old
			if dc_analysis.convergence_check(x=(idiode, vd), dx=(di, dvd), residuum=((vd-ports_v[0])/RS+idiode, ports_v[0]-vd-idiode*self.RS), nv_minus_one=1)[0]:
				break
		dev.last_vd = vd
		return idiode
	def _safe_exp(self, x):
		return math.exp(x) if x<70 else math.exp(70)+10*x
			
	def _get_i(self, v):
		i = self.IS*(self._safe_exp(v/(self.N*self.VT))-1) \
		    + self.ISR*(self._safe_exp(v/(self.NR*self.VT))-1)
		return i
	def get_gm(self, op_index, ports_v, port_index, dev, rs=True):
		if dev.T != self.T:
			self.set_temperature(dev.T)
		gm = self.IS/(self.N*self.VT)*\
			self._safe_exp(ports_v[0]/(self.N*self.VT)) +\
			self.ISR/(self.NR*self.VT)*\
			self._safe_exp(ports_v[0]/(self.NR*self.VT))
		if rs:
			gm = 1./(self.RS + 1./gm)
		return dev.AREA*gm
	def __str__(self):
		pass
	def set_temperature(self, T):
		#print "T:%g => Eg: %g, IS: %g, BV:%g, RS:%g" % (self.T, self.EG, self.IS, self.BV, self.RS)
		T = float(T)
		self.EG = constants.si.Eg(T)
		ni = constants.si.ni(T)
		self.IS = self.IS*(T/self.T)**(self.XTI/self.N) * math.exp(-constants.e\
			*constants.si.Eg(300)/(self.N*constants.k*T)*(1-T/self.T))
		self.BV = self.BV - self.TBV*(T-self.T)
		self.RS = self.RS*(1+self.TRS*(T - self.T))
		self.T = T
		#print "T:%g => Eg: %g, IS: %g, BV:%g, RS:%g" % (self.T, self.EG, self.IS, self.BV, self.RS)
