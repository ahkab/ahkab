# -*- coding: iso-8859-1 -*-
# pz.py
# Pole-Zero evaluation methods
# Copyright 2014 Giuseppe Venturini

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
This module offers the functions needed to perform a pole-zero evaluation.

Ref.

The Generalized Eigenproblem: Pole-Zero Computation

IMPLEMENTATION OF POLE-ZERO ANALYSIS IN SPICE BASED
ON THE MD METHOD

"""

import numpy as np
import circuit
import dc_analysis 
import devices
import transient

def enlarge_matrix(M):
    if M is None:
        return np.mat(np.zeros((1, 1)))
    else:
        return np.mat(np.vstack((np.hstack((M, np.zeros((M.shape[0], 1)))),
                          np.zeros((1, M.shape[1] + 1)))))
def calculate_poles(mc):
    return calculate_singularities(mc, input_source=None, output_port=None, 
                                   calc_zeros=False, MNA=None, shift=0)[0]

def calculate_singularities(mc, input_source=None, output_port=None, calc_zeros=False, MNA=None, shift=0):
    """Calculate poles and zeros.

    *Parameters:*

    mc : circuit instance
    input_source : string or element
        Ignored if zeros are not being calculated.
    output_port : external node (ref. to gnd) or tuple of external nodes
        Ignored if zeros are not being calculated.
    calc_zeros : bool
        Calculate zeros.
    MNA : ndarray, optional
    shift : float
        Shift frequency at which the algorithm should be run.

    *Returns:*

    (poles, zeros) : tuple of array-like
    
    """
    if calc_zeros:
        if type(input_source) != str:
            input_source = input_source.part_id
        if np.isscalar(output_port):
            output_port = (output_port, mc.gnd)
        output_port = map(mc.ext_node_to_int, output_port)
        RIIN = []
        ROUT = []
    if not MNA:
        MNA, N = dc_analysis.generate_mna_and_N(mc)
    D = transient.generate_D(mc, MNA[1:, 1:].shape)
    MNAinv = np.linalg.inv(MNA[1:, 1:] + shift*D[1:, 1:])
    nodes_m1 = len(mc.nodes_dict) - 1
    vde1 = -1
    MC = np.zeros(N.shape)[1:, :]
    TCM = None
    dei_source = 0
    for e1 in mc:
        if circuit.is_elem_voltage_defined(e1):
            vde1 += 1
        if isinstance(e1, devices.Capacitor):
            MC[e1.n1 - 1, 0] += 1. if e1.n1 > 0 else 0.
            MC[e1.n2 - 1, 0] -= 1. if e1.n2 > 0 else 0.
        elif isinstance(e1, devices.Inductor):
            MC[nodes_m1 + vde1] += -1.
        elif calc_zeros and e1.part_id == input_source:
            if isinstance(e1, devices.VSource):
                MC[nodes_m1 + vde1] += -1.
            elif isinstance(e1, devices.ISource):
                MC[e1.n1 - 1, 0] += 1. if e1.n1 > 0 else 0.
                MC[e1.n2 - 1, 0] -= 1. if e1.n2 > 0 else 0.
            else:
                raise Exception("Unknown input source type %s" % input_source)
        else:
            continue
        TV = -1. * MNAinv * MC
        dei_victim = 0
        vde2 = -1
        for e2 in mc:
            if circuit.is_elem_voltage_defined(e2):
                vde2 += 1
            if isinstance(e2, devices.Capacitor):
                v = 0
                if e2.n1:
                    v += TV[e2.n1 - 1, 0]
                if e2.n2:
                    v -= TV[e2.n2 - 1, 0]
            elif isinstance(e2, devices.Inductor):
                v = TV[nodes_m1 + vde2, 0]                
            else:
                continue
            if calc_zeros and e1.part_id == input_source:
                RIIN += [v]
            else:
                if not dei_source:
                    TCM = enlarge_matrix(TCM)
                TCM[dei_victim, dei_source] = v*e1.value
                dei_victim += 1
        if calc_zeros and e1.part_id == input_source:
            ROUTIN = 0
            o1, o2 = output_port
            if o1:
                ROUTIN += TV[o1 - 1, 0]
            if o2:
               ROUTIN -= TV[o2 - 1, 0]
        else:
            dei_source += 1
        # reset, get ready to restart
        MC[:, :] = 0.
    if TCM is not None:
        if np.linalg.det(TCM):
            poles = 1./(2.*np.pi)*(1./np.linalg.eigvals(TCM) + shift)
        else:
            return calculate_singularities(mc, input_source, output_port, 
                                           calc_zeros, MNA=MNA, 
                                           shift=np.abs(np.random.uniform())*1e3)
    else:
        poles = []
    if calc_zeros and TCM is not None:
        # re-loop, get the ROUT elements
        vde1 = -1
        MC = np.zeros(N.shape)[1:, :]
        ROUT = []
        for e1 in mc:
            if circuit.is_elem_voltage_defined(e1):
                vde1 += 1
            if isinstance(e1, devices.Capacitor):
                MC[e1.n1 - 1, 0] += 1. if e1.n1 > 0 else 0.
                MC[e1.n2 - 1, 0] -= 1. if e1.n2 > 0 else 0.
            elif isinstance(e1, devices.Inductor):
                MC[nodes_m1 + vde1, 0] += -1.
            else:
                continue
            TV = -1. * MNAinv * MC
            v = 0
            o1, o2 = output_port
            if o1:
                v += TV[o1 - 1, 0]
            if o2:
               v -= TV[o2 - 1, 0]
            ROUT += [v*e1.value]
            # reset, get ready to restart
            MC[:, :] = 0.
        # Reshape the matrices and evaluate the zero correction.
        RIIN = np.array(RIIN).reshape((-1, 1))
        RIIN = np.tile(RIIN, (1, RIIN.shape[0]))
        ROUT = np.diag(np.atleast_1d(np.array(ROUT)))
        if ROUT.any() and ROUTIN:
            ZTCM = TCM - np.dot(RIIN, ROUT)/ROUTIN
            if np.linalg.det(ZTCM):
                zeros = 1./(2.*np.pi)*(1./np.linalg.eigvals(ZTCM) + shift)
            else:
                return calculate_singularities(mc, input_source, output_port, 
                                           calc_zeros, MNA=MNA, 
                                           shift=np.abs(np.random.uniform())*1e3)
        else:
            zeros = []
    else:
        zeros = []
    return poles, zeros
