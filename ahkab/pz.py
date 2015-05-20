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
This module offers the functions needed to perform a numeric
pole-zero extraction.

Currently, this module implements the MD algorithm, more may be
added in the future.

A description of the algorithm is found in the following references:


    Haley, S.B., "The generalized eigenproblem: pole-zero computation,"
    *Proceedings of the IEEE*, vol.76, no.2, pp.103,120, Feb 1988

and:

    Raghuram, R.; Divekar, D.; Wang, P., "Implementation of pole-zero
    analysis in SPICE based on the MD method," Circuits and Systems, 1991.,
    *Proceedings of the 34th Midwest Symposium on*, pp.380,
    383 vol.1, 14-17 May 1991

Frequency sweeping -- or shifting -- is performed with a random frequency
kick, currently, hoping not to kick so hard that we end up on the negative side.
A bisection method would be better and hopefully will be implemented soon.

Overview
~~~~~~~~

Two main methods are available in this module:

 * :func:`calculate_singularities`, which computes both zeros and poles,
 * :func:`calculate_poles`, which only computes the poles.

Currently this module uses dense matrices.

Reference
~~~~~~~~~

"""
from __future__ import (unicode_literals, absolute_import,
                        division, print_function)

import copy

import numpy as np

from . import circuit
from . import dc_analysis 
from . import devices
from . import transient
from . import plotting
from . import printing
from . import py3compat
from . import options
from . import results

specs = {'pz': {'tokens': ({
                           'label': 'output',
                           'pos': 0,
                           'type': str,
                           'needed': False,
                           'dest': 'output_port',
                           'default': None
                           },
                           {
                           'label': 'input',
                           'pos': 1,
                           'type': str,
                           'needed': False,
                           'dest': 'input_source',
                           'default': None
                           },
                           {
                           'label': 'shift',
                           'pos': 2,
                           'type': float,
                           'needed': False,
                           'dest': 'shift',
                           'default': 0.0
                           })
               }
        }

def _enlarge_matrix(M):
    if M is None:
        return np.zeros((1, 1))
    else:
        return np.vstack((np.hstack((M, np.zeros((M.shape[0], 1)))),
                         np.zeros((1, M.shape[1] + 1))))

def calculate_poles(mc, MNA=None, x0=None, outfile=None, verbose=0):
    """Calculate the circuit poles.

    **Parameters:**

    mc : circuit instance
        The circuit to be analyzed.
    MNA : ndarray, optional
        The Modified Nodal Analysis matrix, if available.
        In case the circuit is non-linear, MNA should include the contributes
        of the non-linear elements (ie the Jacobian :math:`J`). 
    x0 : ndarray or op_solution, optional
         The linearization point. Only needed for non-linear circuits.
    outfile : str or None, optional
        The data filename.
    verbose : int, optional
        Verbosity level, from 0 (silent, default) to 6 (debug).

    **Returns:**

    pz_sol : pz_solution instance
        The PZ solution, with no zeros.
    """
    return calculate_singularities(mc, input_source=None, output_port=None, 
                                   MNA=None, shift=0)[0]

def calculate_singularities(mc, input_source=None, output_port=None, MNA=None,
                            x0=None, shift=0, outfile=None, verbose=0):
    """Calculate poles and zeros.

    By default, only poles are calculated, as they need no information
    other than the circuit description.

    To activate zeros calculation, it is necessary:

    * to specify an input source (``input_source``),
    * to specify an output port (``output_port``).

    **Parameters:**

    mc : circuit instance
        The circuit to be analyzed.
    input_source : string or element, optional
        If zeros are to be calculated, set this to the input surce.
    output_port : external node (ref. to gnd) or tuple of external nodes, opt
        If zeros are to be calculated, set this to the output nodes.
    MNA : ndarray, optional
        The Modified Nodal Analysis matrix, if available.
        In case the circuit is non-linear, MNA should include the contributes
        of the non-linear elements (ie the Jacobian :math:`J`). 
    x0 : ndarray or op_solution, optional
         The linearization point. Only needed for non-linear circuits.
    shift : float, optional
        Shift frequency at which the algorithm should be run.
    outfile : str or None, optional
        The data filename.
    verbose : int, optional
        Verbosity level, from 0 (silent, default) to 6 (debug).

    **Returns:**

    pz_sol : pz_solution instance
        The PZ solution

    """
    calc_zeros = (input_source is not None) and (output_port is not None)
    if calc_zeros:
        if type(input_source) not in py3compat.string_types:
            input_source = input_source.part_id
        if type(output_port) in py3compat.string_types:
            output_port = plotting._split_netlist_label(output_port)[0] 
            output_port = [o[1:] for o in output_port]
            output_port = [o.lower() for o in output_port]
            if len(output_port) == 1:
                # we refer to the ground implicitely
                output_port += ['0']
        if np.isscalar(output_port):
            output_port = (output_port, mc.gnd)
        if (type(output_port) == tuple or type(output_port) == list) \
           and type(output_port[0]) in py3compat.string_types:
            output_port = [mc.ext_node_to_int(o) for o in output_port]
        we_got_source = False
        for e in mc:
            if e.part_id == input_source:
                we_got_source = True
                break
        if not we_got_source:
            raise ValueError('Source %s not found in circuit.' % input_source)
        RIIN = []
        ROUT = []
    if MNA is None:
        MNA, N = dc_analysis.generate_mna_and_N(mc)
        if mc.is_nonlinear():
            # setup x0
            if x0 is None:
                printing.print_warning("PZ: No linearization point provided. Using x0 = 0.")
                x0 = np.zeros((MNA.shape[0] - 1, 1))
            else:
                if isinstance(x0, results.op_solution):
                    x0 = x0.asarray()
                # else
                    # hopefully x0 is an ndarray!
                printing.print_info_line(("Using the supplied op as " +
                                      "linearization point.", 5), verbose)
            J, _ = dc_analysis.build_J_and_Tx(x0, MNA.shape[0]-1, mc, time=0., sparse=False)
            MNA[1:, 1:] += J
    D = transient.generate_D(mc, MNA[1:, 1:].shape)
    MNAinv = np.linalg.inv(MNA[1:, 1:] + shift*D[1:, 1:])
    nodes_m1 = mc.get_nodes_number() - 1
    vde1 = -1
    MC = np.zeros((MNA.shape[0] - 1, 1))
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
        TV = -1. * np.dot(MNAinv, MC)
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
                    TCM = _enlarge_matrix(TCM)
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
            return calculate_singularities(mc, input_source, output_port, MNA=MNA,
                                           outfile=outfile, 
                                           shift=shift+np.abs(1+np.random.uniform())*1e3)
    else:
        poles = []
    if calc_zeros and TCM is not None:
        # re-loop, get the ROUT elements
        vde1 = -1
        MC = np.zeros((MNA.shape[0] - 1, 1))
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
            TV = -1.*np.dot(MNAinv, MC)
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
        if ROUT.any():
            try:
                if not np.all(ROUTIN) or not np.isfinite(ROUTIN):
                    # immediate catch
                    raise ValueError("ROUT-IN is either Inf, NaN or null")
                ZTCM = TCM - np.dot(RIIN, ROUT)/ROUTIN
                if not (np.isfinite(ZTCM).all()):
                    raise ValueError("Array contains infs, NaNs or both")
                    # immediate catch
                eigvals = np.linalg.eigvals(ZTCM)
                if not np.all(eigvals) or not np.isfinite(eigvals).all():
                    # immediate catch
                    raise ValueError("ZTCM eigenvalues contain either Inf, NaN or null values")
                zeros = 1./(2.*np.pi)*(1./np.linalg.eigvals(ZTCM) + shift)
            except ValueError:
                return calculate_singularities(mc, input_source, output_port, 
                                           MNA=MNA, outfile=outfile,
                                           shift=shift+np.abs(np.random.uniform()+1)*1e3)
        elif shift < 1e12:
            return calculate_singularities(mc, input_source, output_port, 
                                       MNA=MNA, outfile=outfile,
                                       shift=shift*np.abs(np.random.uniform()+1)*10)
        else:
            zeros = []
    else:
        zeros = []
    poles = np.array([a for a in poles if np.abs(a) < options.pz_max], dtype=np.complex_)
    zeros = np.array([a for a in zeros if np.abs(a) < options.pz_max], dtype=np.complex_)
    poles = np.sort_complex(poles)
    zeros = np.sort_complex(zeros)
    res = results.pz_solution(mc, poles, zeros, outfile)
    return res

def _check_singularities(res, ref, atol=1e-4, rtol=1e-3):
    ref = copy.copy(ref)
    for i in res:
        success = False
        for si, s in enumerate(ref):
            if np.allclose(i, s, atol=atol, rtol=rtol):
                success = True
                break
            assert si != len(ref) - 1
        ref.pop(si)
    assert not len(ref)

