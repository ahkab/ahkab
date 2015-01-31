# -*- coding: iso-8859-1 -*-
# ac.py
# AC analysis module
# Copyright 2010 Giuseppe Venturini
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

""" This module contains the methods required to perform an AC analysis.

Our AC analysis problem can be written as:

.. math::

    MNA\\ x + AC(\\omega )\\ x + J x + N_{ac}(\\omega) = 0

We need:

    1. the mna matrix :math:`MNA`,
    2. the :math:`AC` matrix, holding the frequency dependent parts,
    3. :math:`J`, the Jacobian matrix from the linearized non-linear elements,
    4. :math:`N_{ac}`, the AC sources contribution.

An OP has to be computed first if there is any non-linear device in the circuit.

When all the matrices are available, it is possible to solve the system
for the frequency values specified by the user, providing the resulting
matrix is not singular (and possibly well conditioned). 

**Building the AC matrix:**

It's easy to set up the voltage lines, since line 2 refers to
node 2, etc...

A capacitor between two example nodes ``n1`` and ``n2`` introduces the
following elements:

.. math::

    \\mathrm{(KCL\\ node\\ n1)}\\qquad +j\\omega C\\ V(n1) - j\\omega C V(n2) + ... = ...

    \\mathrm{(KCL\\ node\\ n2)}\\qquad -j\\omega C\\ V(n1) + j\\omega C V(n2) + ... = ...

Inductors generate, together with voltage sources, ccvs, vcvs, a
additional line in the :math:`MNA` matrix, and hence in :math:`AC` too.
The current flowing through the device gets added to the unknowns vector,
:math:`x`.

For example, in the case of an inductors, we have:

.. math::

    \\mathrm{(KVL\\ over\\ n1\\ and\\ n2)}\\qquad V(n1) - V(n2) - j\\omega L\\ I(\\mathrm{inductor}) = 0

To understand on which line is the KVL line for an inductor, we use the
*order* of the elements in :mod:`ahkab.circuit`:
First are assembled all the voltage rows, then the current ones in the same order in which
the elements that introduce them are found in :mod:`ahkab.circuit`.


Module reference
----------------

"""

from __future__ import (unicode_literals, absolute_import,
                        division, print_function)

import sys

import numpy as np

from . import (circuit, dc_analysis, devices, options, printing, results,
               utilities)

specs = {'ac': {'tokens': ({
                           'label': 'type',
                           'pos': 0,
                           'type': str,
                           'needed': False,
                           'dest': 'sweep_type',
                           'default': options.ac_log_step
                           },
                          {
                           'label': 'start',
                           'pos': 1,
                           'type': float,
                           'needed': True,
                           'dest': 'start',
                           'default': None
                          },
                          {
                           'label': 'stop',
                           'pos': 2,
                           'type': float,
                           'needed': True,
                           'dest': 'stop',
                           'default': None
                          },
                          {
                           'label': 'nsteps',
                           'pos': 3,
                           'type': float,
                           'needed': True,
                           'dest': 'points',
                           'default': None
                          })
               }
        }


def ac_analysis(circ, start, points, stop, sweep_type=None,
                x0=None, mna=None, AC=None, Nac=None, J=None,
                outfile="stdout", verbose=3):
    """Performs an AC analysis.

    **Parameters:**

    circ : Circuit instance
        The circuit to be simulated.

    start : float
        The start angular frequency for the AC analysis

    points : float,
        The number of points to be use the discretize the
        ``[start, stop]`` interval.

    stop : float
        The stop angular frequency.

    sweep_type : string, optional
        Either ``options.ac_log_step`` (``'LOG'``) or ``options.ac_lin_step``
        (``'LIN'``), defaults to a logarithmic sweep.

    x0 : op results instance, optional
        The linearization point. If not set, it will be computed
        running an OP analysis.

    mna, AC, Nax, J : ndarrays, optional
        The matrices to perform the analysis. They will be computed
        if not supplied.

    outfile : string, optional
        the filename of the output file where the results will be written.
        '.ac' is automatically added at the end to prevent different
        analyses from overwriting each-other's results.
        If unset or set to ``None``, defaults to ``stdout``.

    verbose : int, optional
        the verbosity level, from 0 (silent) to 6 (debug).

    **Returns:**

    ACresult : AC solution
        The AC analysis results
    """

    if outfile == 'stdout':
        verbose = 0

    # check step/start/stop parameters
    nsteps = points - 1
    if start == 0:
        printing.print_general_error("AC analysis has start frequency = 0")
        sys.exit(5)
    if start > stop:
        printing.print_general_error("AC analysis has start > stop")
        sys.exit(1)
    if nsteps < 1:
        printing.print_general_error("AC analysis has number of steps <= 1")
        sys.exit(1)
    if sweep_type == options.ac_log_step or sweep_type is None:
        omega_iter = utilities.log_axis_iterator(stop, start, nsteps)
    elif sweep_type == options.ac_lin_step:
        omega_iter = utilities.lin_axis_iterator(stop, start, nsteps)
    else:
        printing.print_general_error("Unknown sweep type.")
        sys.exit(1)

    tmpstr = "Vea =", options.vea, "Ver =", options.ver, "Iea =", options.iea, "Ier =", \
        options.ier, "max_ac_nr_iter =", options.ac_max_nr_iter
    printing.print_info_line((tmpstr, 5), verbose)
    del tmpstr

    printing.print_info_line(("Starting AC analysis: ", 1), verbose)
    tmpstr = "w: start = %g Hz, stop = %g Hz, %d steps" % (start, stop, nsteps)
    printing.print_info_line((tmpstr, 3), verbose)
    del tmpstr

    # It's a good idea to call AC with prebuilt MNA matrix if the circuit is
    # big
    if mna is None:
        (mna, N) = dc_analysis.generate_mna_and_N(circ, verbose=verbose)
        del N
        mna = utilities.remove_row_and_col(mna)
    if Nac is None:
        Nac = _generate_Nac(circ)
        Nac = utilities.remove_row(Nac, rrow=0)
    if AC is None:
        AC = _generate_AC(circ, [mna.shape[0], mna.shape[0]])
        AC = utilities.remove_row_and_col(AC)

    if circ.is_nonlinear():
        if J is not None:
            pass
            # we used the supplied linearization matrix
        else:
            if x0 is None:
                printing.print_info_line(
                    ("Starting OP analysis to get a linearization point...", 3), verbose, print_nl=False)
                # silent OP
                x0 = dc_analysis.op_analysis(circ, verbose=0)
                if x0 is None:  # still! Then op_analysis has failed!
                    printing.print_info_line(("failed.", 3), verbose)
                    printing.print_general_error(
                        "OP analysis failed, no linearization point available. Quitting.")
                    sys.exit(3)
                else:
                    printing.print_info_line(("done.", 3), verbose)
            printing.print_info_line(
                ("Linearization point (xop):", 5), verbose)
            if verbose > 4:
                x0.print_short()
            printing.print_info_line(
                ("Linearizing the circuit...", 5), verbose, print_nl=False)
            J = _generate_J(xop=x0.asmatrix(), circ=circ, mna=mna,
                           Nac=Nac, data_filename=outfile, verbose=verbose)
            printing.print_info_line((" done.", 5), verbose)
            # we have J, continue
    else:  # not circ.is_nonlinear()
        # no J matrix is required.
        J = 0

    printing.print_info_line(("MNA (reduced):", 5), verbose)
    printing.print_info_line((mna, 5), verbose)
    printing.print_info_line(("AC (reduced):", 5), verbose)
    printing.print_info_line((AC, 5), verbose)
    printing.print_info_line(("J (reduced):", 5), verbose)
    printing.print_info_line((J, 5), verbose)
    printing.print_info_line(("Nac (reduced):", 5), verbose)
    printing.print_info_line((Nac, 5), verbose)

    sol = results.ac_solution(circ, start=start, stop=stop,
                              points=nsteps+1, stype=sweep_type, op=x0, 
                              outfile=outfile)

    # setup the initial values to start the iteration:
    nv = len(circ.nodes_dict)
    j = np.complex('j')

    Gmin_matrix = dc_analysis.build_gmin_matrix(
        circ, options.gmin, mna.shape[0], verbose)

    iter_n = 0  # contatore d'iterazione
    printing.print_info_line(("Solving... ", 3), verbose, print_nl=False)

    x = x0
    for omega in omega_iter:
        (x, error, solved, n_iter) = dc_analysis.dc_solve(
            mna=(mna + np.multiply(j * omega, AC) + J),
            Ndc = Nac,
            Ntran = 0,
            circ = circuit.Circuit(
                title="Dummy circuit for AC", filename=None),
            Gmin = Gmin_matrix,
            x0 = x,
            time = None,
            locked_nodes = None,
            MAXIT = options.ac_max_nr_iter,
            skip_Tt = True,
            verbose = 0)
        if solved:
            iter_n = iter_n + 1
            # hooray!
            sol.add_line(omega, x)
        else:
            break

    if solved:
        printing.print_info_line(("done.", 1), verbose)
        ret_value = sol
    else:
        printing.print_info_line(("failed.", 1), verbose)
        ret_value = None

    return ret_value


def _generate_AC(circ, shape):
    """Generates the AC coefficients matrix.

    **Parameters:**

    circ : Circuit instance
        The circuit instance for which the matrix will be generated.

    shape : int
        The reduced MNA size.

    **Returns:**
 
    AC : ndarray
        the *unreduced* AC matrix

    """
    AC = np.matrix(np.zeros((shape[0] + 1, shape[1] + 1)))
    nv = len(circ.nodes_dict)  # - 1
    i_eq = 0  # each time we find a vsource or vcvs or ccvs, we'll add one to this.
    for elem in circ:
        if circuit.is_elem_voltage_defined(elem) and not isinstance(elem, devices.Inductor):
            i_eq = i_eq + 1
        elif isinstance(elem, devices.Capacitor):
            n1 = elem.n1
            n2 = elem.n2
            AC[n1, n1] = AC[n1, n1] + elem.value
            AC[n1, n2] = AC[n1, n2] - elem.value
            AC[n2, n2] = AC[n2, n2] + elem.value
            AC[n2, n1] = AC[n2, n1] - elem.value
        elif isinstance(elem, devices.Inductor):
            AC[nv + i_eq, nv + i_eq] = -1 * elem.value
            if len(elem.coupling_devices):
                for cd in elem.coupling_devices:
                    # get `part_id` of the other inductor (eg. "L32")
                    other_id_wdescr = cd.get_other_inductor(elem.part_id)
                    # find its index to know which column corresponds to its
                    # current
                    other_index = circ.find_vde_index(
                        other_id_wdescr, verbose=0)
                    # add the term.
                    AC[nv + i_eq, nv + other_index] += -1 * cd.M
            i_eq = i_eq + 1

    if options.cmin > 0:
        cmin_mat = np.matrix(np.eye(shape[0] + 1 - i_eq))
        cmin_mat[0, 1:] = 1
        cmin_mat[1:, 0] = 1
        cmin_mat[0, 0] = cmin_mat.shape[0] - 1
        if i_eq:
            AC[:-i_eq, :-i_eq] += options.cmin * cmin_mat
        else:
            AC += options.cmin * cmin_mat

    return AC


def _generate_Nac(circ):
    """Generate the vector holding the contribution of AC sources.

    **Parameters:**

    circ : Circuit instance
        The circuit instance for which the matrix will be generated.
    """
    n_of_nodes = len(circ.nodes_dict)
    Nac = np.mat(np.zeros((n_of_nodes, 1)), dtype=complex)
    j = np.complex('j')
    # process `ISource`s
    for elem in circ:
        if isinstance(elem, devices.ISource) and elem.abs_ac is not None:
            # convenzione normale!
            Nac[elem.n1, 0] = Nac[elem.n1, 0] + \
                elem.abs_ac * np.exp(j * elem.arg_ac)
            Nac[elem.n2, 0] = Nac[elem.n2, 0] - \
                elem.abs_ac * np.exp(j * elem.arg_ac)
    # process vsources
    # for each vsource, introduce a new variable: the current flowing through it.
    # then we introduce a KVL equation to be able to solve the circuit
    for elem in circ:
        if circuit.is_elem_voltage_defined(elem):
            index = Nac.shape[0]
            Nac = utilities.expand_matrix(Nac, add_a_row=True, add_a_col=False)
            if isinstance(elem, devices.VSource) and elem.abs_ac is not None:
                Nac[index, 0] = -1.0 * elem.abs_ac * np.exp(j * elem.arg_ac)
    return Nac


def _generate_J(xop, circ, mna, Nac, data_filename, verbose=0):
    """Build the linearized matrix :math:`J`.
    """
    # setup J
    J = np.mat(np.zeros(mna.shape))
    Tlin = np.mat(np.zeros(Nac.shape))
    for elem in circ:
        if elem.is_nonlinear:
            dc_analysis.update_J_and_Tx(J, Tlin, xop, elem, time=None)
    # del Tlin # not needed! **DC**!

    return J
