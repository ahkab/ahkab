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

- mosq_device
- mosq_model

The Square Law Mos Model

Assuming, :math:`V_{ds} > 0` and a transistor type N in the
following, we have the following regions implemented:

1. No subthreshold conduction.
   :math:`V_{gs} < V_T`
   :math:`I_D = 0`
2. Ohmic region
   :math:`V_{GS} > V_T` and :math:`V_{GD} > V_T`
   :math:`I_D = k_n W/L ((V_{GS}-V_{T})V_{DS} - V_{DS}^2/2)`
3. Saturation region
   :math:`V_{GS} > V_T` and :math:`V_{DS} > V_{GS} - V_{T}`
   :math:`V_{GS} < V_{T}`
   :math:`I_D = 1/2 k_n W/L (V_{GS}-V_T)^2 * [1 + \lambda*(V_{DS}-V_{GS}+V_T)]`

"""
from __future__ import division

from __future__ import (unicode_literals, absolute_import,
                        division, print_function)

import math
import numpy as np

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
LAMBDA_DEFAULT = .5
AVT_DEFAULT = 7.1e-3 * 1e-6
AKP_DEFAULT = 1.8e-2 * 1e-6

TCV_DEFAULT = 1e-3
BEX_DEFAULT = -1.5

ISMALL_GUESS_MIN = 1e-10


class mosq_device:
    def __init__(self, part_id, nd, ng, ns, nb, W, L, model, M=1, N=1):
        """Quadratic Law MOSFET device

        **Parameters:**

        part_id : string
            The part ID of the model. Eg. ``'M1'`` or ``'Mlow'``, the first
            letter should always be ``'M'``.
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
        model : mosq_mos_model instance
            the model for the device
        M : int, optional
            shunt multiplier (n. of shunt devices)
        N : int, optional
            series multiplier (n. of series devices)

        Selected methods:
        - get_output_ports() -> (nd, ns)
        - get_drive_ports() -> (nd, nb), (ng, nb), (ns, nb)
        """
        self.ng = ng
        self.nb = nb
        self.n1 = nd
        self.n2 = ns
        self.ports = ((self.n1, self.n2), (
            self.ng, self.n2), (self.nb, self.n2))

        class dev_class:
            pass  # empty class to hold device parameters
        self.device = dev_class()
        self.device.L = float(L)  # channel length -
        self.device.W = float(W)  # channel width -
        self.device.M = int(M)  # parallel multiple device number
        self.device.N = int(N)  # series multiple device number
        self.device.mckey = None
        self.device.part_id = part_id
        self.mosq_model = model
        self.mc_enabled = False
        self.opdict = {}
        self.opdict.update(
            {'state': (float('nan'), float('nan'), float('nan'))})
        self.part_id = part_id
        self.is_nonlinear = True
        self.is_symbolic = True
        self.dc_guess = [self.mosq_model.VTO*0.4*self.mosq_model.NPMOS,
                         self.mosq_model.VTO*1.1*self.mosq_model.NPMOS,
                         0]

        devcheck, reason = self.mosq_model._device_check(self.device)
        if not devcheck:
            raise ValueError(reason + " out of boundaries.")

    def get_drive_ports(self, op):
        """Get the circuit ports that drive the device.

        **Returns:**

        a tuple of tuples of ports nodes, as:

        ``(port0, port1, port2...)``

        Where each port is in the form:

            port0 = (nplus, nminus)

        """
        return self.ports  # d,g,b

    def get_output_ports(self):
        """Get the circuit ports where the device injects current.

        **Returns:**

        ports : a tuple of tuples of ports nodes, as:

        ``(port0, port1, port2...)``

        Where each port is in the form:

            port0 = (nplus, nminus)

        """
        return ((self.n1, self.n2),)

    def __str__(self):
        mos_type = self._get_mos_type()
        rep = self.part_id + " %(nd)s %(ng)s %(ns)s %(nb)%s " + \
              self.mosq_model.name + " w=" + str(self.device.W) + " l=" + \
              str(self.device.L) + " M=" + str(self.device.M) + " N=" + \
              str(self.device.N)
        return rep

    def _get_mos_type(self):
        """Returns N or P (capitalized), depending on the device type.
        """
        mtype = 'N' if self.mosq_model.NPMOS == 1 else 'P'
        return mtype

    def istamp(self, ports_v, time=0, reduced=True):
        """Get the current matrix

        A matrix corresponding to the current flowing in the element
        with the voltages applied as specified in the ``ports_v`` vector.

        **Parameters:**

        ports_v : list
            A list in the form: [voltage_across_port0, voltage_across_port1, ...]
        time: float
            the simulation time at which the evaluation is performed.
            It has no effect here. Set it to None during DC analysis.

        """
        sw_vect, CS = self.mosq_model.get_voltages(*ports_v)
        ids = self.mosq_model.get_ids(self.device, sw_vect, self.opdict)
        istamp = np.array((CS*ids, -CS*ids), dtype=np.float64)
        indices = ((self.n1 - 1*reduced, self.n2 - 1*reduced), (0, 0))
        if reduced:
            delete_i = [pos for pos, i in enumerate(indices[0]) if i == -1]
            istamp = np.delete(istamp, delete_i, axis=0)
            indices = tuple(zip(*[(i, j) for i, j in zip(*indices) if i != -1]))
        return indices, istamp

    def update_status_dictionary(self, ports_v):
        if self.opdict is None:
            self.opdict = {}
        if not (self.opdict['state'] == ports_v[0] and 'gmd' in self.opdict) or \
            not (self.opdict['state'] == ports_v[0] and 'gm' in self.opdict) or \
            not (self.opdict['state'] == ports_v[0] and 'gmb' in self.opdict) or \
                not (self.opdict['state'] == ports_v[0] and 'Ids' in self.opdict):

            swapped, CS = self.mosq_model.get_voltages(*ports_v[0])
            self.opdict['state'] == ports_v[0]
            gstamp = self.gstamp(ports_v[0], reduced=False)[1]
            self.opdict['gmd'] = gstamp[0, 0]
            self.opdict['gm'] = gstamp[0, 1]
            self.opdict['gmb'] = gstamp[0, 3]
            self.opdict['Ids'] = self.istamp(ports_v[0], reduced=False)[1][0]

    def print_op_info(self, ports_v):
        arr = self.get_op_info(ports_v)
        print(arr, end=' ')

    def get_op_info(self, ports_v):
        """Operating point info, for design/verification. """
        self.update_status_dictionary(ports_v)
        sat_status = "SATURATION" if self.opdict['SAT'] else "LINEAR"
        if not self.opdict["ON"]:
            status = "OFF"
        else:
            status = "ON"

        arr = [
            [self.part_id + " ch", status, "", "", sat_status, "", "", "", "", "", "", ""], ]
        arr.append(
            ["beta", "[A/V^2]:", self.opdict['beta'], "Weff", "[m]:", str(self.opdict['W']) + " (" + str(self.device.W) + ")",
             "L", "[m]:", str(self.opdict['L']) + " (" + str(self.device.L) + ")", "M/N:", "", str(self.device.M) + "/" + str(self.device.N)])
        arr.append(["Vds", "[V]:", float(ports_v[0][0]), "Vgs", "[V]:", float(
            ports_v[0][1]), "Vbs", "[V]:", float(ports_v[0][2]),  "", "", ""])
        arr.append(["VTH", "[V]:", self.opdict['VTH'], "VOD", "[V]:", self.opdict[
                   'VOD'], "", "", "", "VA", "[V]:", str(self.opdict['Ids'] / self.opdict['gmd'])])
        arr.append(
            ["Ids", "[A]:", self.opdict['Ids'], "", "", "", "", "", "", "", "", ''])
        arr.append(["gm", "[S]:", self.opdict['gm'], "gmb", "[S]:",
                   self.opdict['gmb'], "ro", "[Ohm]:", 1 / self.opdict['gmd'], "", "", ""])

        return printing.table_setup(arr)

    def gstamp(self, ports_v, time=0, reduced=True):
        """Returns the differential (trans)conductance rs the port specified by port_index
        when the element has the voltages specified in ports_v across its ports,
        at (simulation) time.

        ports_v: a list in the form: [voltage_across_port0, voltage_across_port1, ...]
        port_index: an integer, 0 <= port_index < len(self.get_ports())
        time: the simulation time at which the evaluation is performed. Set it to
        None during DC analysis.
        """
        indices = ([self.n1 - 1]*4 + [self.ng - 1]*4 + [self.n2 - 1]*4 + [self.nb - 1]*4,
                   [self.n1 - 1, self.ng - 1, self.n2 - 1, self.nb - 1]*4)

        sw_vect, CS = self.mosq_model.get_voltages(*ports_v)
        gmd = self.mosq_model.get_gmd(self.device, sw_vect, self.opdict)
        gmg = self.mosq_model.get_gm(self.device, sw_vect, self.opdict)
        gmb = self.mosq_model.get_gmb(self.device, sw_vect, self.opdict)
        if gmd == 0:
            gmd = options.gmin*2
        if gmg == 0:
            gmg = options.gmin*2
        if gmb == 0:
            gmb = -2*options.gmin
        stamp = np.array(((gmd, gmg, -gmd-gmb-gmg, gmb),
                          (0, 0, 0, 0),
                          (-gmd, -gmg, gmd + gmg + gmb, -gmb),
                          (0, 0, 0, 0)), dtype=np.float64)
        if CS == -1:
            stamp = self.mosq_model.T1*stamp*self.mosq_model.T2
        if (self.opdict['state'] != ports_v[0]).any():
            self.opdict = {'state':ports_v[0]}
        self.opdict.update({'gmd': stamp[0, 0]})
        self.opdict.update({'gm': stamp[0, 1]})
        self.opdict.update({'gmb': stamp[0, 3]})
        if reduced:
            zap_rc = [pos for pos, i in enumerate(indices[1][:4]) if i == -1]
            stamp = np.delete(stamp, zap_rc, axis=0)
            stamp = np.delete(stamp, zap_rc, axis=1)
            indices = tuple(zip(*[(i, y) for i, y in zip(*indices) if (i != -1 and y != -1)]))
            stamp_flat = stamp.reshape(-1)
            stamp_folded = []
            indices_folded = []
            for ix, it in enumerate([(i, y) for i, y in zip(*indices)]):
                if not it in indices_folded:
                    indices_folded.append(it)
                    stamp_folded.append(stamp_flat[ix])
                else:
                    w = indices_folded.index(it)
                    stamp_folded[w] += stamp_flat[ix]
            indices = tuple(zip(*indices_folded))
            stamp = np.array(stamp_folded)
        return indices, stamp

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

    def get_netlist_elem_line(self, nodes_dict):
        mos_type = self._get_mos_type()
        return "%s %s %s %s %s %s type=%s w=%g l=%g m=%g n=%g" % \
              (self.part_id, nodes_dict[self.n1], nodes_dict[self.ng],
              nodes_dict[self.n2], nodes_dict[self.nb], self.mosq_model.name,
              mos_type, self.device.W, self.device.L, self.device.M,
              self.device.N)


class mosq_mos_model:

    def __init__(self, name=None, TYPE='n', TNOM=None, COX=None,
                 GAMMA=None, NSUB=None, PHI=None, VTO=None, KP=None,
                 LAMBDA=None, AKP=None, AVT=None,
                 TOX=None, VFB=None, U0=None, TCV=None, BEX=None):

        self.name = "model_mosq0" if name is None else name
        Vth = constants.Vth()
        self.TNOM = float(TNOM) if TNOM is not None else constants.Tref
        # print "TYPE IS:" + TYPE
        if TYPE.lower() == 'n':
            self.NPMOS = 1
        elif TYPE.lower() == 'p':
            self.NPMOS = -1
        else:
            raise ValueError("Unknown MOS type %s" % TYPE)

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
                2 * constants.e * constants.si.esi * NSUB * 10 ** 6 / self.COX)
        else:
            self.GAMMA = GAMMA_DEFAULT
        if PHI is not None:
            self.PHI = float(PHI)
        elif NSUB is not None:
            self.PHI = 2 * constants.Vth(self.TNOM) * math.log(
                NSUB * 10 ** 6 / constants.si.ni(self.TNOM))
        else:
            self.PHI = PHI_DEFAULT
        if VTO is not None:
            self.VTO = self.NPMOS * float(VTO)
            if self.VTO < 0:
                print("(W): model %s has internal negative VTO (%f V)." % (self.name, self.VTO))
        elif VFB is not None:
            self.VTO = VFB + PHI + GAMMA * PHI  # inv here??
        else:
            self.VTO = VTO_DEFAULT

        if KP is not None:
            self.KP = float(KP)
        elif U0 is not None:
            self.KP = (U0 * 10 ** -4) * self.COX
        else:
            self.KP = KP_DEFAULT

        self.LAMBDA = LAMBDA if LAMBDA is not None else LAMBDA_DEFAULT
        # Intrinsic model temperature parameters
        self.TCV = self.NPMOS * \
            float(TCV) if TCV is not None else self.NPMOS * TCV_DEFAULT
        self.BEX = float(BEX) if BEX is not None else BEX_DEFAULT

        # Monte carlo
        self.AVT = AVT if AVT is not None else AVT_DEFAULT
        self.AKP = AKP if AKP is not None else AKP_DEFAULT

        self.set_device_temperature(constants.T)

        sc, sc_reason = self._self_check()
        if not sc:
            raise Exception(sc_reason + " out of range")
        self.T1 = np.array(((0, 0, 1, 0),
                            (0, 1, 0, 0),
                            (1, 0, 0, 0),
                            (0, 0, 0, 1)))
        self.T2 = np.array(((0, 0, 1, 0),
                            (0, 1, 0, 0),
                            (1, 0, 0, 0),
                            (0, 0, 0, 1)))

    def set_device_temperature(self, T):
        """Change the temperature of the device. VTO, KP and PHI get updated.
        """
        self.TEMP = T
        self.VTO = self.VTO - self.TCV * (T - self.TNOM)
        self.KP = self.KP * (T / self.TNOM) ** self.BEX
        self.PHI = self.PHI * T / self.TNOM + 3 * constants.Vth(self.TNOM) * math.log(T / self.TNOM) \
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
            [self.name, "", "", TYPE + " MOS", "SQUARE MODEL", "", "", "", "",  "", "", ""])
        arr.append(["KP", "[A/V^2]", self.KP, "VTO", "[V]:", self.VTO,
                   "TOX", "[m]", self.TOX, "COX", "[F/m^2]:", self.COX])
        arr.append(["PHI", "[V]:", self.PHI, "GAMMA", "sqrt(V)", self.GAMMA,
                   "NSUB", "[cm^-3]", self.NSUB,  "VFB", "[V]:", self.VFB])
        arr.append(
            ["U0", "[cm^2/(V*s)]:", self.U0, "TCV", "[V/K]", self.TCV, "BEX", "", self.BEX,  "", "", ""])
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
        vds = vds * self.NPMOS
        vgs = vgs * self.NPMOS
        vbs = vbs * self.NPMOS
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
        # print ((float(vds_new), float(vgs_new), float(vbs_new)), cs)
        return (float(vds_new), float(vgs_new), float(vbs_new)), cs

    def get_svt_skp(self, device, debug=False):
        if device.mckey and debug:
            print("Monte carlo enabled. key:", device.mckey)
        if device.mckey:
            svt = device.mckey[0] * self.AVT / math.sqrt(
                2 * device.W * device.L)
            skp = device.mckey[1] * self.AKP / math.sqrt(
                2 * device.W * device.L)
        else:
            svt, skp = 0, 0
        return svt, skp

    def get_ids(self, device, xxx_todo_changeme, opdict=None, debug=False):
        """Returns:
            IDS, the drain-to-source current (de-normalized),
            qs, the (scaled) charge at the source,
            qr, the (scaled) charge at the drain.
        """
        (vds, vgs, vbs) = xxx_todo_changeme
        if debug:
            print("=== %s (%sch) current for vds: %g, vgs: %g, vbs: %g" \
                  % (device.part_id, 'n'*(self.NPMOS == 1) +
                  'p'*(self.NPMOS == -1), vds, vgs, vbs))

        # monte carlo support
        svt, skp = self.get_svt_skp(device, debug=debug)

        if debug:
            print("PHI:", self.PHI, "vbs:", vbs)

        vsqrt1 = max(-vbs + 2*self.PHI, 0.)
        vsqrt2 = max(2*self.PHI, 0.)

        VT = self.VTO + svt + self.GAMMA * \
            (math.sqrt(vsqrt1) - math.sqrt(vsqrt2))
        if vgs < VT:
            ids = options.iea * (vgs / VT + vds / VT) / 100
            if debug:
                print("OFF: %g" % ids)
        else:
            if vds < vgs - VT -0.5*self.LAMBDA*(VT - vgs)**2:
                ids = (skp + 1) * self.KP * device.W / \
                      device.L * ((vgs - VT) * vds - .5 * vds ** 2)
                if debug:
                    print("OHMIC: %g" % ids)
            else:
                ids = (skp + 1) * .5 * self.KP * device.W / device.L * (
                      vgs - VT) ** 2 * (1 + self.LAMBDA * (vds - vgs + VT + 0.25*self.LAMBDA*(VT - vgs)**2))
                if debug:
                    print("SAT: %g" % ids)
        Ids = self.NPMOS * device.M / device.N * ids

        opdict.update(
            {'state': (vds * self.NPMOS, vgs * self.NPMOS, vbs * self.NPMOS)})
        opdict.update(
            {'Ids': Ids, "W": device.W, "L": device.L, "ON": 1 * (vgs >= VT)})
        opdict.update({'beta': .5 * self.KP * device.W / device.L})
        opdict.update(
            {'VTH': VT, "VOD": self.NPMOS * (vgs - VT), 'SAT': vds > vgs - VT})

        return Ids

        if debug:
            print("vd:", vd, "vg:", VG / self.scaling.Ut, "vs:", vs, "vds:", vd - vs)
            print("v_ifn:", v_ifn, "v_irn:", v_irn)
            print("ifn:", ifn, "irn:", irn)
            print("ip_abs_err:", ip_abs_err)
            print("Vth:", self.scaling.Ut)
            print("nv", nv, "Is", self.scaling.Is)
            print("Weff:", device.W, "Leff:", Leff)
            print("NPMOS:", self.NPMOS, "CS_FACTOR", CS_FACTOR)

    def get_gmb(self, device, xxx_todo_changeme1, opdict=None, debug=False):
        """Returns the source-bulk transconductance or d(IDS)/d(VS-VB)."""
        (vds, vgs, vbs) = xxx_todo_changeme1
        svt, skp = self.get_svt_skp(device, debug=False)
        assert vds >= 0
        vsqrt1 = max(-vbs + 2*self.PHI, 0.)
        vsqrt2 = max(2*self.PHI, 0.)
        VT = self.VTO + svt + self.GAMMA * \
            (math.sqrt(vsqrt1) - math.sqrt(vsqrt2))
        gmb = 0
        if vgs < VT:
            pass # gmb = 0
        else:
            if vds < vgs - VT:
                if vsqrt1 > 0:
                    gmb = self.KP * self.GAMMA * vds * device.W / \
                          (2 * device.L * vsqrt1 ** .5)
            else:
                if vsqrt1 > 0:
                    gmb += -0.25*self.KP*self.GAMMA*self.LAMBDA*device.W * (vsqrt1 > 0) * \
                           (-self.GAMMA*(-vsqrt2**.5 + vsqrt1**.5) + vgs - self.VTO)**2 / \
                           (device.L * vsqrt1**.5)
                    gmb += +0.5*self.KP*self.GAMMA*device.W*(self.LAMBDA* \
                            (self.GAMMA * (vsqrt2**.5 + vsqrt1**.5) + vds - vgs + self.VTO) + 1.0) *\
                            (-self.GAMMA * (vsqrt2**.5 + vsqrt1**.5) \
                            + vgs - self.VTO) / (device.L * vsqrt1**.5)
        gmb = self.NPMOS * (1 + skp) * gmb * device.M / device.N
        if debug:
            print("gmb %g" % gmb)
        return gmb

    def get_gmd(self, device, xxx_todo_changeme2, opdict=None, debug=False):
        """Returns the drain-bulk transconductance or d(IDS)/d(VD-VB)."""
        (vds, vgs, vbs) = xxx_todo_changeme2
        svt, skp = self.get_svt_skp(device, debug=False)
        assert vds >= 0
        vsqrt1 = max(-vbs + 2*self.PHI, 0.)
        vsqrt2 = max(2*self.PHI, 0.)
        VT = self.VTO + svt + self.GAMMA * \
            (math.sqrt(vsqrt1) - math.sqrt(vsqrt2))
        if vgs < VT:
            gmd = options.iea / VT / 100
        else:
            if vds < vgs -VT -0.5*self.LAMBDA*(VT - vgs)**2: # correction term disc. due to LAMBDA
                gmd = self.KP * device.W / device.L * (vgs - vds - VT)
            else:
                gmd = 0.5 * self.KP * self.LAMBDA * device.W / device.L * \
                      (vgs - VT)**2
        gmd = (1 + skp) * gmd * device.M / device.N
        if debug:
            print("gmd %g" % gmd)
        return gmd

    def get_gm(self, device, xxx_todo_changeme3, opdict=None, debug=False):
        """Returns the gate-bulk transconductance or d(IDS)/d(VG-VB)."""
        (vds, vgs, vbs) = xxx_todo_changeme3
        svt, skp = self.get_svt_skp(device, debug=False)
        assert vds >= 0
        vsqrt1 = max(-vbs + 2*self.PHI, 0.)
        vsqrt2 = max(2*self.PHI, 0.)
        VT = self.VTO + svt + self.GAMMA * \
            (math.sqrt(vsqrt1) - math.sqrt(vsqrt2))
        if vgs < VT:
            gm = options.iea / VT / 100
        else:
            if vds < vgs - VT:
                gm = self.KP * device.W / device.L * vds
            else:
                gm = -0.5*self.KP*self.LAMBDA * device.W/device.L * (-self.GAMMA*(-vsqrt2**.5 + vsqrt1**.5) + vgs - self.VTO)**2 \
                     +0.5*self.KP * device.W/device.L *(self.LAMBDA*( self.GAMMA*(-vsqrt2**.5 + vsqrt1**.5) + vds - vgs + self.VTO) + 1.0) *\
                    (-2 * self.GAMMA * (-vsqrt2**.5 + vsqrt1**.5) + 2*vgs - 2*self.VTO)
        gm = (1 + skp) * gm * device.M / device.N
        if debug:
            print("gmg %g" % gm)
        return gm

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
        elif self.AVT and self.AVT < 0:
            ret = (False, "AVT " + str(self.AVT))
        elif self.AKP and self.AKP < 0:
            ret = (False, "AKP " + str(self.AKP))
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

