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
Partial implementation of the EKV 3.0 MOS transistor model

The EKV model was developed by Matthias Bucher, Christophe Lallement,
Christian Enz, Fabien Théodoloz, François Krummenacher at the Electronics
Laboratories, Swiss Federal Institute of Technology (EPFL),
Lausanne, Switzerland.

The Tecnical Report upon which this implementation
is based is available here:

`EKV Technical Report <http://legwww.epfl.ch/ekv/pdf/ekv_v262.pdf>`_.

This module defines two classes:

* :class:`ekv_device`
* :class:`ekv_mos_model`


**Features:**

* EKV model implementation, computation of charges, potentials,
  reverse and forward currents, slope factor and normalization factors.
* Calculation of trans-conductances based on the charge-driven approach.
* N/P MOS symmetry
* Rudimentary temperature effects.

**The Missing Features:**

* Channel length modulation,
* Reverse Short Channel Effect (RSCE),
* Complex mobility degradation,
* Transcapacitances,
* Quasi-static implementation,

Patches to implement the above are welcome!

.. note::
    The default values in the model are suitable for a generic
    500nm feature-size CMOS process.

"""

from __future__ import (unicode_literals, absolute_import,
                        division, print_function)

import scipy, scipy.optimize
import math

from scipy.optimize import newton

from . import constants
from . import options
from . import utilities
from . import printing


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
    """EKV device

    **Parameters:**

    part_id : string
        The element identifier, eg 'M1'

    nd : int
        drain node
    ng : int
        gate node
    ns : int
        source node
    nb : int
        bulk node
    L : float
        element width [m]
    W : float
        element length [m]
    M : int
        multiplier (n. of shunt devices)
    N : int
        series mult. (n. of series devices)
    model : ekv_model instance
        The corresponding instance of ekv_mos_model

    Selected methods:
    - get_output_ports() -> (nd, ns)
    - get_drive_ports() -> (nd, nb), (ng, nb), (ns, nb)

    """

    INIT_IFRN_GUESS = 1

    def __init__(self, part_id, nd, ng, ns, nb, W, L, model, M=1, N=1):
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
            raise Exception(reason + " out of boundaries.")

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
        if not (self.opdict['state'] == ports_v[0] and 'gmd' in self.opdict) or \
            not (self.opdict['state'] == ports_v[0] and 'gmg' in self.opdict) or \
            not (self.opdict['state'] == ports_v[0] and 'gms' in self.opdict) or \
                not (self.opdict['state'] == ports_v[0] and 'Ids' in self.opdict):

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

    def get_op_info(self, ports_v):
        """Information regarding the Operating Point (OP)

        **Parameters:**

        ports_v : list of lists
            The voltages applied to all the driving ports, grouped by output
            port.

        i.e.

        ::

            [<list of voltages for the drive ports of output port 0>,
             <list of voltages for the drive ports of output port 1>,
             ...,
             <list of voltages for the drive ports of output port N>]

        Usually, this method returns ``op_keys`` and the corresponding
        ``op_info``, two lists, one holding the labels, the other the
        corresponding values.

        In the case of MOSFETs, the values are way too many to be shown in a
        linear table. For this reason, we return ``None`` as ``op_keys``, and we
        return for ``op_info`` a list which holds both labels and values in a
        table-like manner, spanning the vertical and horizontal dimension.

        For this reason, each MOSFET has to have its OP info printed alone, not
        grouped as it happens with most other elements.

        **Returns:**

        op_keys : ``None``
            See above for why this value is always ``None``.
        op_info : list of floats
            The OP information ready to be passed to :func:`printing.table` for
            arranging it in a pretty table to display.
        """
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
                   self.opdict['gms'], "rob", u"[\u2126]:", 1 / self.opdict['gmd'], "", "", ""])
        arr.append(
            ["if:", "", self.opdict['ifn'], "ir:", "", self.opdict['irn'],
             "Qf", "[C/m^2]:", self.opdict["qf"], "Qr", "[C/m^2]:", self.opdict["qr"], ])
        # arr.append([  "", "", "", "", "", ""])

        return None, arr

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

    def get_netlist_elem_line(self, nodes_dict):
        mos_type = self._get_mos_type()
        return "%s %s %s %s %s %s type=%s w=%g l=%g m=%g n=%g" % \
              (self.part_id, nodes_dict[self.n1], nodes_dict[self.ng],
              nodes_dict[self.n2], nodes_dict[self.nb], self.ekv_model.name,
              mos_type, self.device.W, self.device.L, self.device.M,
              self.device.N)


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
            self.GAMMA = math.sqrt(2.0*constants.e*constants.si.esi*NSUB
                                   *10**6/self.COX)
        else:
            self.GAMMA = GAMMA_DEFAULT
        if PHI is not None:
            self.PHI = float(PHI)
        elif NSUB is not None:
            self.PHI = 2. * constants.Vth(self.TNOM) * \
                       math.log(NSUB*10.0**6.0/constants.si.ni(self.TNOM))
        else:
            self.PHI = PHI_DEFAULT
        if VTO is not None:
            self.VTO = self.NPMOS * float(VTO)
            if self.VTO < 0:
                printing.print_warning("model %s has internal negative VTO (%f V)."
                                       % (self.name, self.VTO))
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
            raise Exception(sc_reason + " out of range")

    def set_device_temperature(self, T):
        """Change the temperature of the device.

        Correspondingly, ``VTO``, ``KP`` and ``PHI`` get updated.
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
        print(printing.table(arr))

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

    def get_ids(self, device, xxx_todo_changeme, opdict=None, debug=False):
        """Returns:
            IDS, the drain-to-source current (de-normalized),
            qs, the (scaled) charge at the source,
            qr, the (scaled) charge at the drain.
        """
        (vd, vg, vs) = xxx_todo_changeme
        if debug:
            print("=== Current for vd:", vd, "vg:", vg, "vs:", vs)
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
            print("Scaled voltages: vd:", vd, "vp:", vp, "vs:", vs)

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
            print("vd:", vd, "vg:", VG / self.scaling.Ut, "vs:", vs, "vds:", vd - vs)
            print("v_ifn:", v_ifn, "v_irn:", v_irn)
            print("ifn:", ifn, "irn:", irn)
            print("ip_abs_err:", ip_abs_err)
            print("Vth:", self.scaling.Ut)
            print("nv", nv, "Is", self.scaling.Is)
            print("Weff:", device.W, "Leff:", Leff)
            print("NPMOS:", self.NPMOS, "CS_FACTOR", CS_FACTOR)

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
            print("current:", Ids)

        return Ids, qf, qr

    def get_leq_virp(self, device, xxx_todo_changeme1, Vp, Leff, ifn):
        # if ifn > 0 and Vp - constants.Vth()*vd > 0:
        (vd, vg, vs) = xxx_todo_changeme1
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

    def get_gms(self, device, xxx_todo_changeme2, opdict=None, debug=False):
        """Returns the source-bulk transconductance or d(IDS)/d(VS-VB)."""
        (vd, vg, vs) = xxx_todo_changeme2
        (j1, j2, j3), CS_FACTOR = self.get_voltages(vd, vg, vs)
        Ids, qf, qr = self.get_ids(device, (vd, vg, vs), opdict, debug)
        if CS_FACTOR == +1:
            gms = -1.0 * self.scaling.Gs * qf
        elif CS_FACTOR == -1:
            gms = -self.scaling.Gs * qr
        return gms

    def get_gmd(self, device, xxx_todo_changeme3, opdict=None, debug=False):
        """Returns the drain-bulk transconductance or d(IDS)/d(VD-VB)."""
        (vd, vg, vs) = xxx_todo_changeme3
        (j1, j2, j3), CS_FACTOR = self.get_voltages(vd, vg, vs)
        Ids, qf, qr = self.get_ids(device, (vd, vg, vs), opdict, debug)
        if CS_FACTOR == +1:
            gmd = self.scaling.Gs * qr
        elif CS_FACTOR == -1:
            gmd = self.scaling.Gs * qf
        return gmd

    def get_gmg(self, device, xxx_todo_changeme4, opdict=None, debug=False):
        """Returns the gate-bulk transconductance or d(IDS)/d(VG-VB)."""
        (vd, vg, vs) = xxx_todo_changeme4
        VP, nv, nq = self.get_vp_nv_nq(float(vg))
        Ids, qf, qr = self.get_ids(device, (vd, vg, vs), opdict, debug)
        (j1, j2, j3), CS_FACTOR = self.get_voltages(vd, vg, vs)
        gmg = CS_FACTOR * self.scaling.Gs * (qf - qr) / nv
        return gmg

    def get_ismall(self, vsmall, ip_abs_err, iguess=None, debug=False):
        """Solves the problem: given v, find i such that:

        .. math::
            v = ln(q) + 2q

        ..math::
            q = sqrt(.25 + i) - .5

        The Newton Method is used inside.
        """
        # starting guess for Newton's Method.
        if iguess is None:
            iguess = 1.0
        # sanity checks
        if math.isnan(vsmall):
            raise ValueError("Attempted to calculate a current corresponding to a NaN voltage.")
        if not ip_abs_err > 0:
            raise ValueError("The normalized current absolute error has been set to a negative value.")
        ismall = newton(self._vsmall_obj, iguess, fprime=self._vsmall_obj_prime,
                        fprime2=self._vsmall_obj_prime2, args=(vsmall,),
                        tol=1.48e-8, maxiter=500)
        #print(ismall, max(ismall, 0))
        return max(ismall, 10*utilities.EPS)

    def _vsmall_obj(self, x, vsmall):
        """Returns :math:`e` according to the equations:
            q = sqrt(.25 + x) - .5
            e = ln(q) + 2q - vsmall
        """
        if x > 0:
            return math.log(math.sqrt(.25 + x) - 0.5) + 2.0 * math.sqrt(.25 + x) - 1.0 - vsmall
        else:
            return x - vsmall

    def _vsmall_obj_prime(self, x, vsmall):
        """The Newton algorithm in get_ismall(...) requires the evaluation of the
        first derivative of the fixed point function:
            dv/di = 1.0/(sqrt(.25+i)-.5) * .5/sqrt(.25 + i) + 1/sqrt(.25 + i)
        This is provided by this module.
        """
        if x < utilities.EPS:
            x = 10*utilities.EPS
        return 0.5/(math.sqrt(.25 + x) - .5)/math.sqrt(.25 + x) + \
               1.0/math.sqrt(.25 + x)

    def _vsmall_obj_prime2(self, x, vsmall):
        if x < utilities.EPS:
            x = 10*utilities.EPS
        return -1./(4*(x + 0.25)*(math.sqrt(x + 0.25) - 0.5)**2) - 1./(2.*(x + 0.25)**(3./2.)) - 1./(4.*(x + 0.25)**(3./2.)*(math.sqrt(x + 0.25) - 0.5))

    def get_vsmall(self, ismall, verbose=3):
        """Returns v according to the equations:
            q = sqrt(.25 + i) - .5
            v = ln(q) + 2q
        """
        if abs(ismall) < utilities.EPS:
            ismall = utilities.EPS  # otherwise we get log(0)
            if verbose == 6:
                print("EKV: Machine precision limited the resolution on i. (i<EPS)")
        vsmall = math.log(math.sqrt(.25 + ismall) - 0.5) + \
                 2.0 * math.sqrt(.25 + ismall) - 1.0
        return vsmall

    def get_dvsmall_dismall(self, ismall, verbose=3):
        """The Newton algorithm in get_ismall(...) requires the evaluation of the
        first derivative of the fixed point function:
            dv/di = 1.0/(sqrt(.25+i)-.5) * .5/sqrt(.25 + i) + 1/sqrt(.25 + i)
        This is provided by this module.
        """
        if abs(ismall) < utilities.EPS:
            ismall = utilities.EPS
            if verbose == 6:
                print("EKV: Machine precision limited the resolution on dv/di in the NR iteration. (i<EPS)")
        dvdi = 1.0 / (math.sqrt(.25 + ismall) - .5) * .5 / \
            math.sqrt(.25 + ismall) + 1.0 / math.sqrt(.25 + ismall)
        return dvdi

    def ismall2qsmall(self, ismall, verbose=0):
        """ i(f,r) -> q(f,r)
        Convert a source/drain scaled current to the corresponding normalized charge."""
        if verbose == 6:  # ismall is lower than EPS, errors here are usually not important
            print("EKV: Machine precision limited the resolution on q(s,d). (i<EPS)")
        ismall = max(0, ismall)
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

