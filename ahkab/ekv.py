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
- Calculation of trans-conductances based on the charge approach.
- N/P MOS symmetry
- Rudimentary temperature effects.

The Missing Features:
- Channel length modulation
- Reverse Short Channel Effect (RSCE)
- Complex mobility degradation is missing
- Transcapacitances
- Quasistatic implementation
"""

__version__ = "0.091"

import constants
import options
import utilities
import printing
import math


# DEFAULT VALUES FOR 500n CH LENGTH
COX_DEFAULT = .7e-3
VTO_DEFAULT = .5
GAMMA_DEFAULT = 1.0
PHI_DEFAULT = .7
KP_DEFAULT = 50e-6
UCRIT_DEFAULT = 2e6
LAMBDA_DEFAULT = .5
XJ_DEFAULT = .1e-6


TCV_DEFAULT = 1e-3
BEX_DEFAULT = -1.5

ISMALL_GUESS_MIN = 1e-10


class ekv_device:
    INIT_IFRN_GUESS = 1

    def __init__(self, nd, ng, ns, nb, W, L, model, M=1, N=1, part_id='M'):
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
        self.ports = ((self.n1, self.nb), (
            self.ng, self.nb), (self.n2, self.nb))

        class dev_class:
            pass  # empty class to hold device parameters
        self.device = dev_class()
        self.device.L = float(L)  # channel length -
        self.device.W = float(W)  # channel width -
        self.device.M = int(M)  # parallel multiple device number
        self.device.N = int(N)  # series multiple device number
        self.ekv_model = model
        self.opdict = {}
        self.opdict.update(
            {'state': (float('nan'), float('nan'), float('nan'))})
        self.opdict.update({'ifn': self.INIT_IFRN_GUESS})
        self.opdict.update({'irn': self.INIT_IFRN_GUESS})
        self.opdict.update(
            {'ip_abs_err': self.ekv_model.get_ip_abs_err(self.device)})
        self.part_id = part_id
        self.is_nonlinear = True
        self.is_symbolic = True
        self.dc_guess = [self.ekv_model.VTO * (
            0.1) * self.ekv_model.NPMOS, self.ekv_model.VTO * (1.1) * self.ekv_model.NPMOS, 0]

        devcheck, reason = self.ekv_model._device_check(self.device)
        if not devcheck:
            raise Exception, reason + " out of boundaries."

    def get_drive_ports(self, op):
        """Returns a tuple of tuples of ports nodes, as:
        (port0, port1, port2...)
        Where each port is in the form:
        port0 = (nplus, nminus)
        """
        return self.ports  # d,g,s

    def get_output_ports(self):
        return ((self.n1, self.n2),)

    def __str__(self):
        mos_type = self._get_mos_type()
        rep = " " + self.ekv_model.name + " w=" + str(self.device.W) + " l=" + \
            str(self.device.L) + " M=" + str(self.device.M) + " N=" + \
            str(self.device.N)

        return rep

    def _get_mos_type(self):
        """Returns N or P (capitalized)
        """
        mtype = 'N' if self.ekv_model.NPMOS == 1 else 'P'
        return mtype

    def i(self, op_index, ports_v, time=0):
        """Returns the current flowing in the element with the voltages
        applied as specified in the ports_v vector.

        ports_v: [voltage_across_port0, voltage_across_port1, ...]
        time: the simulation time at which the evaluation is performed.
              It has no effect here. Set it to None during DC analysis.

        """
        ret, j1, j2 = self.ekv_model.get_ids(self.device, ports_v,
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
            TEF = abs(gms * constants.Vth() / ids)
        self.opdict['TEF'] = TEF

    def print_op_info(self, ports_v):
        arr = self.get_op_info(ports_v)
        print arr,

    def get_op_info(self, ports_v):
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

        arr = [
            [self.part_id, mos_type.upper() + " ch", wmsi_status, "", "", sat_status, "", "", "", "", "", ""], ]
        arr.append(
            ["beta", "[A/V^2]:", self.opdict['beta'], "Weff", "[m]:", str(self.opdict['Weff']) + " (" + str(self.device.W) + ")",
             "Leff", "[m]:", str(self.opdict['Leff']) + " (" + str(self.device.L) + ")", "M/N:", "", str(self.device.M) + "/" + str(self.device.N)])
        arr.append(
            ["Vdb", "[V]:", float(ports_v[0][0]), "Vgb", "[V]:", float(ports_v[0][1]),
             "Vsb", "[V]:", float(ports_v[0][2]),  "Vp", "[V]:", self.opdict['Vp'], ])
        arr.append(
            ["VTH", "[V]:", self.opdict['VTH'], "VOD", "[V]:", self.opdict['VOD'],
             "nq: ", "", self.opdict['nq'], "VA", "[V]:", str(self.opdict['Ids'] / self.opdict['gmd'])])
        arr.append(
            ["Ids", "[A]:", self.opdict['Ids'], "nv: ", "", self.opdict['nv'],
             "Ispec", "[A]:", self.opdict["Ispec"], "TEF:", "", str(self.opdict['TEF']), ])
        arr.append(["gmg", "[S]:", self.opdict['gmg'], "gms", "[S]:",
                   self.opdict['gms'], "rob", "[Ohm]:", 1 / self.opdict['gmd'], "", "", ""])
        arr.append(
            ["if:", "", self.opdict['ifn'], "ir:", "", self.opdict['irn'],
             "Qf", "[C/m^2]:", self.opdict["qf"], "Qr", "[C/m^2]:", self.opdict["qr"], ])
        # arr.append([  "", "", "", "", "", ""])

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
            g = sign * options.gmin * 2.0

        # print type(g), g

        if op_index == 0 and port_index == 0:
            self.opdict.update({'gmd': g})
        elif op_index == 0 and port_index == 1:
            self.opdict.update({'gmg': g})
        elif op_index == 0 and port_index == 2:
            self.opdict.update({'gms': g})

        return g

    def get_value_function(self, identifier):
        def get_value(self):
            return self.opdict[identifier]
        return get_value


class scaling_holder:
    pass  # will hold the scaling factors


class ekv_mos_model:

    def __init__(self, name=None, TYPE='n', TNOM=None, COX=None,
                 GAMMA=None, NSUB=None, PHI=None, VTO=None, KP=None,
                 XJ=None, LAMBDA=None,
                 TOX=None, VFB=None, U0=None, TCV=None, BEX=None):

        self.scaling = scaling_holder()

        self.name = "model_ekv0" if name is None else name
        self.TNOM = float(TNOM) if TNOM is not None else constants.Tref
        self.NPMOS = 1 if TYPE == 'n' else -1

        # optional parameters (no defaults)
        self.TOX = float(TOX) if TOX is not None else None
        self.NSUB = float(NSUB) if NSUB is not None else None
        self.VFB = self.NPMOS * float(VFB) if VFB is not None else None
        self.U0 = float(U0) if U0 is not None else None

        # crucial parameters
        if COX is not None:
            self.COX = float(COX)
        elif TOX is not None:
            self.COX = constants.si.eox / TOX
        else:
            self.COX = COX_DEFAULT

        if GAMMA is not None:
            self.GAMMA = float(GAMMA)
        elif NSUB is not None:
            self.GAMMA = math.sqrt(
                2.0 * constants.e * constants.si.esi * NSUB * 10 ** 6 / self.COX)
        else:
            self.GAMMA = GAMMA_DEFAULT
        if PHI is not None:
            self.PHI = float(PHI)
        elif NSUB is not None:
            self.PHI = 2. * constants.Vth(self.TNOM) * math.log(
                NSUB * 10.0 ** 6.0 / constants.si.ni(self.TNOM))
        else:
            self.PHI = PHI_DEFAULT
        if VTO is not None:
            self.VTO = self.NPMOS * float(VTO)
            if self.VTO < 0:
                print "(W): model %s has internal negative VTO (%f V)." % (self.name, self.VTO)
        elif VFB is not None:
            self.VTO = VFB + PHI + GAMMA * PHI  # inv here??
        else:
            self.VTO = self.NPMOS * VTO_DEFAULT

        if KP is not None:
            self.KP = float(KP)
        elif U0 is not None:
            self.KP = (U0 * 10.0 ** -4) * self.COX
        else:
            self.KP = KP_DEFAULT

        self.LAMBDA = LAMBDA if LAMBDA is not None else LAMBDA_DEFAULT
        self.XJ = XJ if XJ is not None else XJ_DEFAULT
        self.UCRIT = UCRIT_DEFAULT
        # Intrinsic model temperature parameters
        self.TCV = self.NPMOS * \
            float(TCV) if TCV is not None else self.NPMOS * TCV_DEFAULT
        self.BEX = float(BEX) if BEX is not None else BEX_DEFAULT

        self.set_device_temperature(constants.T)

        # Setup switches
        self.SATLIM = math.exp(4.0)
        self.WMSI_factor = 10
        self.NR_damp_factor = options.nl_voltages_lock_factor

        sc, sc_reason = self._self_check()
        if not sc:
            raise Exception, sc_reason + " out of range"

    def set_device_temperature(self, T):
        """Change the temperature of the device. VTO, KP and PHI get updated.
        """
        self.TEMP = T
        self.VTO = self.VTO - self.TCV * (T - self.TNOM)
        self.KP = self.KP * (T / self.TNOM) ** self.BEX
        self.PHI = self.PHI * T / self.TNOM + 3.0 * constants.Vth(self.TNOM) * math.log(T / self.TNOM) \
            - constants.si.Eg(
                self.TNOM) * T / self.TNOM + constants.si.Eg(T)

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
        arr.append(
            [self.name, "", "", TYPE + " MOS", "EKV MODEL", "", "", "", "",  "", "", ""])
        arr.append(["KP", "[A/V^2]", self.KP, "VTO", "[V]:", self.VTO,
                   "TOX", "[m]", self.TOX, "COX", "[F/m^2]:", self.COX])
        arr.append(["PHI", "[V]:", self.PHI, "GAMMA", "sqrt(V)", self.GAMMA,
                   "NSUB", "[cm^-3]", self.NSUB,  "VFB", "[V]:", self.VFB])
        arr.append(
            ["U0", "[cm^2/(V*s)]:", self.U0, "TCV", "[V/K]", self.TCV, "BEX", "", self.BEX,  "", "", ""])
        arr.append(["INTERNAL", "", "", "SAT LIMIT", "", self.SATLIM,
                   "W/M/S INV FACTOR", "", self.WMSI_factor,  "", "", ""])
        printing.table_print(arr)

    def get_voltages(self, vd, vg, vs):
        """Performs the VD <-> VS swap if needed.
        Returns:
        (VD, VG, VS) after the swap
        CS, an integer which equals to:
            +1 if no swap was necessary,
            -1 if VD and VS have been swapped.
        """
        # vd / vs swap
        vd = vd * self.NPMOS
        vg = vg * self.NPMOS
        vs = vs * self.NPMOS
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
        """Absolute error to be enforced in the calculation of the normalized currents.
        """
        return options.iea / (2.0 * constants.Vth(self.TEMP) ** 2.0 * self.KP * device.M * device.W / device.L)

    def setup_scaling(self, nq, device):
        """Calculates and stores in self.scaling the following factors:
          Ut, the thermal voltage,
          Is, the specific current,
          Gs, the specific transconductance,
          Qs, the specific charge.
        """
        self.scaling.Ut = constants.Vth()
        self.scaling.Is = 2 * nq * \
            self.scaling.Ut ** 2 * self.KP * device.W / device.L
        self.scaling.Gs = 2 * nq * \
            self.scaling.Ut * self.KP * device.W / device.L
        self.scaling.Qs = 2 * nq * self.scaling.Ut * self.COX
        return

    def get_vp_nv_nq(self, VG):
        """Calculates and returns:
            VP, the pinch-off voltage,
            nv, the slope factor,
            nq, the charge linearization factor.
        """
        VGeff = VG - self.VTO + self.PHI + self.GAMMA * math.sqrt(self.PHI)
        if VGeff > 0 and VG - self.VTO + (math.sqrt(self.PHI) + self.GAMMA / 2.0) ** 2 > 0:
            VP = VG - self.VTO - self.GAMMA * \
                (math.sqrt(VG - self.VTO + (math.sqrt(self.PHI) + self.GAMMA / 2.0) ** 2)
                 - (math.sqrt(self.PHI) + self.GAMMA / 2.0))
            if math.isnan(VP):
                VP = 0  # the argument of sqrt ^^ went negative
        else:
            VP = -self.PHI
        # print "VG", VG, "VGeff", VGeff, "VP", VP, self.GAMMA, self.PHI,
        # math.sqrt(VG -self.VTO +(math.sqrt(self.PHI)+self.GAMMA/2)**2), VG
        # -self.VTO +(math.sqrt(self.PHI)+self.GAMMA/2)**2
        nq = 1.0 + .5 * self.GAMMA / math.sqrt(self.PHI + .5 * VP)
        nv = 1.0 + .5 * self.GAMMA / math.sqrt(self.PHI + VP + 1e-12)

        return VP, nv, nq

    def get_ids(self, device, (vd, vg, vs), opdict=None, debug=False):
        """Returns:
            IDS, the drain-to-source current (de-normalized),
            qs, the (scaled) charge at the source,
            qr, the (scaled) charge at the drain.
        """
        if debug:
            print "=== Current for vd:", vd, "vg:", vg, "vs:", vs
        ip_abs_err = self.get_ip_abs_err(device) if opdict[
            'ip_abs_err'] is None else opdict['ip_abs_err']

        (VD, VG, VS), CS_FACTOR = self.get_voltages(vd, vg, vs)

        # Weff, Leff = self.get_eff_wl(device.W, device.L)

        VP, nv, nq = self.get_vp_nv_nq(VG)

        self.setup_scaling(nq, device)

        vp = VP / self.scaling.Ut
        vs = VS / self.scaling.Ut
        vd = VD / self.scaling.Ut

        if debug:
            print "Scaled voltages: vd:", vd, "vp:", vp, "vs:", vs

        v_ifn = vp - vs
        ifn = self.get_ismall(v_ifn, opdict['ip_abs_err'], max(
            opdict['ifn'], ISMALL_GUESS_MIN), debug=debug)

        if False:
            Leff = device.L
            v_irn = vp - vd
        else:
            Leff, v_irn = self.get_leq_virp(
                device, (vd, vg, vs), VP, device.L, ifn)

        irn = self.get_ismall(v_irn, opdict['ip_abs_err'], max(
            opdict['irn'], ISMALL_GUESS_MIN), debug=debug)

        if debug:
            print "vd:", vd, "vg:", VG / self.scaling.Ut, "vs:", vs, "vds:", vd - vs
            print "v_ifn:", v_ifn, "v_irn:", v_irn
            print "ifn:", ifn, "irn:", irn
            print "ip_abs_err:", ip_abs_err
            print "Vth:", self.scaling.Ut
            print "nv", nv, "Is", self.scaling.Is
            print "Weff:", device.W, "Leff:", Leff
            print "NPMOS:", self.NPMOS, "CS_FACTOR", CS_FACTOR

        qf = self.ismall2qsmall(ifn)
        qr = self.ismall2qsmall(irn)

        Ids =  CS_FACTOR * self.NPMOS * device.L / \
            Leff * device.M * self.scaling.Is * (ifn - irn)

        vd_real = vd if CS_FACTOR == 1 else vs
        vs_real = vs if CS_FACTOR == 1 else vd

        opdict.update(
            {'state': (vd_real * self.NPMOS, vg * self.NPMOS, vs_real * self.NPMOS)})
        opdict.update({'Ids': Ids, "Weff": device.W, "Leff": Leff, 'Vp': VP})
        opdict.update({'ifn': ifn, "irn": irn, "nv": nv, "nq": nq,
                      'beta': .5 * self.KP * device.W / Leff, 'Ispec': self.scaling.Is})
        opdict.update(
            {'VTH': self.VTO, "VOD": self.NPMOS * nv * (VP - VS), 'SAT': ifn > irn * self.SATLIM})
        opdict.update({'qf': qf * self.scaling.Qs, 'qr': qr * self.scaling.Qs})

        if max(ifn, irn) > self.WMSI_factor:
            WMSI = 2
        elif max(ifn, irn) < 1. / self.WMSI_factor:
            WMSI = 0
        else:
            WMSI = 1
        opdict.update({'WMSI': WMSI})

        if debug:
            print "current:", Ids

        return Ids, qf, qr

    def get_leq_virp(self, device, (vd, vg, vs), Vp, Leff, ifn):
        # if ifn > 0 and Vp - constants.Vth()*vd > 0:
        assert vd >= vs
        Vc = self.UCRIT * device.N * Leff
        Vdss  = Vc * \
            (math.sqrt(.25 + constants.Vth() / Vc * math.sqrt(ifn)) - .5)
             # eq. 46
        # Drain-to-source saturation voltage for reverse normalized current,
        # eq. 47
        Vdssp = Vc * (math.sqrt(.25 + constants.Vth() / Vc * (math.sqrt(ifn) - .75 * math.log(ifn))) - .5) + \
            constants.Vth() * (math.log(.5 * Vc / constants.Vth()) - .6)

        # channel length modulation
        vser_1 = math.sqrt(ifn) - Vdss / constants.Vth()
        # if vser_1 < 0:
        #   vser_1 = 0
        Vds = (vd - vs) * .5 * constants.Vth()
        delta_v = 4.0 * constants.Vth() * math.sqrt(
            self.LAMBDA * vser_1 + 1.0 / 64)  # eq. 48
        Vip = math.sqrt(Vdss ** 2 + delta_v ** 2) - math.sqrt(
            (Vds - Vdss) ** 2 + delta_v ** 2)  # eq 50
        Lc = math.sqrt(constants.si.esi * self.XJ / self.COX)  # eq. 51
        delta_l = self.LAMBDA * Lc * \
            math.log(1 + (Vds - Vip) / (Lc * self.UCRIT))  # eq. 52

        # Equivalent channel length including channel-length modulation and
        # velocity saturation
        Lp = device.N * Leff - delta_l + (Vds + Vip) / self.UCRIT  # eq. 53
        Lmin = device.N * Leff / 10.0  # eq. 54
        Leq = .5 * (Lp + math.sqrt(Lp ** 2 + Lmin ** 2))  # eq. 55

        assert not math.isnan(Vdssp)
        assert not math.isnan(delta_v)

        v_irp = (Vp - Vds - vs * constants.Vth() - math.sqrt(Vdssp ** 2 + delta_v ** 2)
                 + math.sqrt((Vds - Vdssp) ** 2 + delta_v ** 2)) / constants.Vth()
        # else:
        #   v_irp = Vp/constants.Vth() - vd
        #   Leq = Leff

        return Leq, v_irp

    def get_gms(self, device, (vd, vg, vs), opdict=None, debug=False):
        """Returns the source-bulk transconductance or d(IDS)/d(VS-VB)."""
        (j1, j2, j3), CS_FACTOR = self.get_voltages(vd, vg, vs)
        Ids, qf, qr = self.get_ids(device, (vd, vg, vs), opdict, debug)
        if CS_FACTOR == +1:
            gms = -1.0 * self.scaling.Gs * qf
        elif CS_FACTOR == -1:
            gms = -self.scaling.Gs * qr
        return gms

    def get_gmd(self, device, (vd, vg, vs), opdict=None, debug=False):
        """Returns the drain-bulk transconductance or d(IDS)/d(VD-VB)."""
        (j1, j2, j3), CS_FACTOR = self.get_voltages(vd, vg, vs)
        Ids, qf, qr = self.get_ids(device, (vd, vg, vs), opdict, debug)
        if CS_FACTOR == +1:
            gmd = self.scaling.Gs * qr
        elif CS_FACTOR == -1:
            gmd = self.scaling.Gs * qf
        return gmd

    def get_gmg(self, device, (vd, vg, vs), opdict=None, debug=False):
        """Returns the gate-bulk transconductance or d(IDS)/d(VG-VB)."""
        VP, nv, nq = self.get_vp_nv_nq(float(vg))
        Ids, qf, qr = self.get_ids(device, (vd, vg, vs), opdict, debug)
        (j1, j2, j3), CS_FACTOR = self.get_voltages(vd, vg, vs)
        gmg = CS_FACTOR * self.scaling.Gs * (qf - qr) / nv
        return gmg

    def get_ismall(self, vsmall, ip_abs_err, iguess=None, debug=False):
        """Solves the problem: given v, find i such that:
            v = ln(q) + 2q
            q = sqrt(.25 + i) - .5
        A damped Newton algorithm is used inside.
        """
        # starting guess for Newton's Method.
        if iguess is None:
            iguess = 1.0
        # sanity checks
        if math.isnan(vsmall):
            raise Exception, \
                "Attempted to calculate a current corresponding to a NaN voltage."
        if not ip_abs_err > 0:
            raise Exception, \
                "The normalized current absolute error has been set to a negative value."
        # if vsmall < 0:
        #   return 0.0

        check = False
        ismall = iguess
        if debug:
            iter_c = 0

        while True:
            if debug:
                iter_c = iter_c + 1

            vsmall_iter, numeric_problem_v = self.get_vsmall(ismall)
            dvdi, numeric_problem_i = self.get_dvsmall_dismall(ismall)
            deltai = (vsmall - vsmall_iter) / dvdi

            numeric_problem = numeric_problem_i or numeric_problem_v
            if debug:
                print "Numeric problem:", numeric_problem
                print "ABS: deltai < ip_abs_err", deltai, "<", ip_abs_err, ":", abs(deltai) < ip_abs_err
                print "REL: deltai < ismall*options.ier", deltai, "<", ismall * options.ier, abs(deltai) < ismall * options.ier
                print deltai, ismall
            # absolute and relative value convergence checks.
            if ((abs(deltai) < ip_abs_err or numeric_problem) and abs(deltai) < ismall * options.ier) or \
                    (abs(deltai) < ip_abs_err * 1e-6 or numeric_problem):
                # To make the algorithm more robust,
                # the convergence check has to be passed twice in a row
                # to reach convergence.
                if not check:
                    check = True
                else:
                    break
            else:
                check = False
            # convergence was not reached, update ismall
            if math.isnan(ismall):
                print "Ismall is NaN!!"
                exit()
            if ismall == 0:
                # this is a sign we went below the machine resolution
                # it makes no sense to iterate there as quantization errors
                # prevent reaching a meaningful result.
                break
            else:
                # Damped Newton with domain restriction: ismall >= 0.
                ratio = deltai / ismall
                if ratio > self.NR_damp_factor:
                    # Do not allow a change in ismall bigger than self.NR_damp_factor
                    # in a single iteration
                    ismall = self.NR_damp_factor * ismall
                elif ratio <= -1:
                    # this would give a negative ismall
                    ismall = 0.1 * ismall
                else:
                    ismall = ismall + deltai
        if debug:
            print str(iter_c) + " iterations."
        return ismall

    def get_vsmall(self, ismall, verbose=3):
        """Returns v according to the equations:
            q = sqrt(.25 + i) - .5
            v = ln(q) + 2q
        """
        if abs(ismall) < utilities.EPS:
            ismall = utilities.EPS  # otherwise we get log(0)
            if verbose == 6:
                print "EKV: Machine precision limited the resolution on i. (i<EPS)"
            numeric_problem = True
        else:
                    numeric_problem = False
        vsmall = math.log(math.sqrt(.25 + ismall) - 0.5) + \
            2.0 * math.sqrt(.25 + ismall) - 1.0
        return vsmall, numeric_problem

    def get_dvsmall_dismall(self, ismall, verbose=3):
        """The Newton algorithm in get_ismall(...) requires the evaluation of the
        first derivative of the fixed point function:
            dv/di = 1.0/(sqrt(.25+i)-.5) * .5/sqrt(.25 + i) + 1/sqrt(.25 + i)
        This is provided by this module.
        """
        if abs(ismall) < utilities.EPS:
            ismall = utilities.EPS
            numeric_problem = True
            if verbose == 6:
                print "EKV: Machine precision limited the resolution on dv/di in the NR iteration. (i<EPS)"
        else:
            numeric_problem = False
        dvdi = 1.0 / (math.sqrt(.25 + ismall) - .5) * .5 / \
            math.sqrt(.25 + ismall) + 1.0 / math.sqrt(.25 + ismall)
        return dvdi, numeric_problem

    def ismall2qsmall(self, ismall, verbose=0):
        """ i(f,r) -> q(f,r)
        Convert a source/drain scaled current to the corresponding normalized charge."""
        if verbose == 6:  # ismall is lower than EPS, errors here are usually not important
            print "EKV: Machine precision limited the resolution on q(s,d). (i<EPS)"
        qsmall = math.sqrt(.25 + ismall) - .5
        return qsmall

    def qsmall2ismall(self, qsmall):
        """ q(f,r) -> i(f,r)
        Convert a source/drain scaled charge to the corresponding normalized current."""
        ismall = qsmall ** 2 + qsmall
        return ismall

    def _self_check(self):
        """Performs sanity check on the model parameters."""
        ret = True, ""
        if self.NSUB is not None and self.NSUB < 0:
            ret = (False, "NSUB " + str(self.NSUB))
        elif self.U0 is not None and not self.U0 > 0:
            ret = (False, "UO " + str(self.U0))
        elif not self.GAMMA > 0:
            ret = (False, "GAMMA " + str(self.GAMMA))
        elif not self.PHI > 0.1:
            ret = (False, "PHI " + str(self.PHI))
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

    ekv_m = ekv_mos_model(TYPE='n', KP=50e-6, VTO=.4)
    ma = ekv_device(1, 2, 3, 4, W=10e-6, L=1e-6, model=ekv_m)
    ma.part_id = 'M1'

    # OP test
    vd = 0.0
    vg = 1.0
    vs = 0.0
    ma.print_op_info(((vd, vg, vs),))
    ekv_m.print_model()

    # gmUt/Ids test
    import mosq
    msq = mosq.mosq(1, 2, 3, kp=50e-6, w=10e-6,
                    l=1e-6, vt=.4, lambd=0.0, mos_type='n')
    data0 = []
    data1 = []
    data2 = []
    data3 = []
    vd = 2.5
    if True:
        vs = 1
        for Vhel in range(1, 2500):
            print ".",
            vg = Vhel / 1000.0
            ma.update_status_dictionary(((vd, vg, 0),))
            data0.append(ma.opdict['Ids'])
            # print "Current for vd", vd, "vg", vg, "vs", vs
            data1.append(ma.opdict['TEF'])
            isq = msq.i((vd, vg, vs),)
            gmsq = msq.g((vd, vg, vs), 0)
            if isq > 0:
                data2.append(gmsq / isq * constants.Vth())
            else:
                data2.append(float('nan'))
            data3.append(isq)
    plt.semilogx(data0, data1, data3, data2)
    plt.title('Transconductance efficiency factor')
    plt.legend(['(GM*UT)/ID'])
    plt.show()
