# -*- coding: iso-8859-1 -*-
# bfpss.py
# Brute-force PSS analysis module
# Copyright 2006 Giuseppe Venturini

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

"""Brute-force periodic steady state analysis module"""

from __future__ import (unicode_literals, absolute_import,
                        division, print_function)

import sys
import numpy as np
import numpy.linalg
import scipy
import scipy.sparse

from . import circuit
from . import dc_analysis
from . import devices
from . import implicit_euler
from . import options
from . import printing
from . import results
from . import ticker
from . import transient
from . import utilities


def bfpss(circ, period, step=None, points=None, autonomous=False, x0=None,
          mna=None, Tf=None, D=None, outfile='stdout',
          vector_norm=lambda v: max(abs(v)), verbose=3):
    """Performs a PSS analysis employing the 'brute-force' algorithm

    The time step is constant and IE will be used as DF.

    **Parameters:**

    circ : Circuit instance
        the circuit to be simulated

    period : float
        the period of the solution

    step : float, optional
        the time step between consecutive points, it will be calculated
        if not provided

    points : int, optional
        the number of points to be used to sample one period, it will be
        calculated if not provided

    autonomous : bool, optional
        Is the circuit clocked or autonomously oscillating?
        With the current implementation, setting ``autonomous=True``
        will result in an exception being raised, autonomous circuits are
        not supported

    x0 : ndarray, optional
        The initial guess to be used. (Experimental, needs work.)

    mna, D, Tf : ndarrays, optional
        The matrices describing the circuit may be supplied to speed up
        the solution, if available. If not supplied, they will be
        automatically calculated.

    vector_norm : function, optional
       The norm to be employed in the convergence checks.
       Defaults to the Inf norm.

    outfile : str, optional
        the output filename. Defaults to ``'stdout'``.

    verbose : int, optional
        Verbosity level on a scale from 0 (silent) to 6 (very verbose).
        The ``verbose`` flag is automatically set is to zero
        if ``datafilename == 'stdout'``

    .. note::

        ``step`` and ``points`` are mutually exclusive options:

        * if ``step`` is specified, the number of points will be
          automatically determined.
        * if ``points`` is set, the step will be automatically determined.
        * if none of them is set, ``options.bfpss_default_points`` will
          be used as value for ``points`` and ``step`` computed accordingly.

    **Returns:**

    sol : results.pss_solution
        The simulation results

    """
    if outfile == "stdout":
        verbose = 0

    printing.print_info_line(
        ("Starting periodic steady state analysis:", 3), verbose)
    printing.print_info_line(("Method: brute-force", 3), verbose)

    if mna is None or Tf is None:
        (mna, Tf) = dc_analysis.generate_mna_and_N(circ, verbose=verbose)
        mna = utilities.remove_row_and_col(mna)
        Tf = utilities.remove_row(Tf, rrow=0)
    elif not mna.shape[0] == Tf.shape[0]:
        printing.print_general_error(
            "mna matrix and N vector have different number of rows.")
        sys.exit(0)

    if D is None:
        D = transient.generate_D(circ, [mna.shape[0], mna.shape[0]])
        D = utilities.remove_row_and_col(D)
    elif not mna.shape == D.shape:
        printing.print_general_error(
            "mna matrix and D matrix have different sizes.")
        raise ValueError

    (points, step) = utilities.check_step_and_points(step, points, period,
                                           options.bfpss_default_points)

    n_of_var = mna.shape[0]
    locked_nodes = circ.get_locked_nodes()
    tick = ticker.ticker(increments_for_step=1)
    sparse = n_of_var*points > options.dense_matrix_limit

    CMAT = _build_CMAT(mna, D, step, points, tick, n_of_var=n_of_var,
                      sparse=sparse, verbose=verbose)

    x = _build_x(mna, step, points, tick, x0=x0, n_of_var=n_of_var,
                 verbose=verbose)

    Tf = _build_Tf(Tf, points, tick, n_of_var=n_of_var, verbose=verbose)

    # time variable component: Tt this is always the same in each iter.
    # So we build it once for all
    # It holds all time-dependent sources (both V/I).
    Tt = _build_Tt(circ, points, step, tick, n_of_var=n_of_var,
                   verbose=verbose)

    # Indices to differentiate between currents and voltages in the
    # convergence check
    nv_indices = []
    ni_indices = []
    nv_1 = circ.get_nodes_number() - 1
    ni = n_of_var - nv_1
    for i in range(points):
        nv_indices += (i * mna.shape[0] * np.ones(nv_1) + \
                      np.arange(nv_1)).tolist()
        ni_indices += (i * mna.shape[0] * np.ones(ni) + \
                      np.arange(nv_1, n_of_var)).tolist()

    converged = False

    printing.print_info_line(("Solving... ", 3), verbose, print_nl=False)
    tick.reset()
    tick.display(verbose > 2)
    if sparse:
        J = scipy.sparse.lil_matrix(CMAT.shape)
    else:
        J = np.zeros(CMAT.shape)
    T = np.zeros((CMAT.shape[0], 1))
    # td is a np array that will hold the damping factors
    td = np.zeros((points, 1))
    iteration = 0  # newton iteration counter

    while True:
        if iteration:  # the first time are already all zeros
            if sparse:
                J = scipy.sparse.lil_matrix(CMAT.shape)
            else:
                J[:, :] = 0
            T[:, 0] = 0
            td[:, 0] = 0
        for index in range(1, points):
            for elem in circ:
                # build all dT(xn)/dxn (stored in J) and T(x)
                if elem.is_nonlinear:
                    oports = elem.get_output_ports()
                    for opindex in range(len(oports)):
                        dports = elem.get_drive_ports(opindex)
                        v_ports = []
                        for dpindex in range(len(dports)):
                            dn1, dn2 = dports[dpindex]
                            # build v: remember we trashed the
                            # 0 row and 0 col of mna -> -1
                            v = 0
                            if dn1:
                                v = v + x[index * n_of_var + dn1 - 1, 0]
                            if dn2:
                                v = v - x[index * n_of_var + dn2 - 1, 0]
                            v_ports.append(v)
                        # all drive ports are ready.
                        n1, n2 = oports[opindex][0], oports[opindex][1]
                        if n1:
                            T[index * n_of_var + n1 - 1, 0] = T[index * n_of_var + n1 - 1, 0] + \
                                                              elem.i(opindex, v_ports)
                        if n2:
                            T[index * n_of_var + n2 - 1, 0] = T[index * n_of_var + n2 - 1, 0] - \
                                                              elem.i(opindex, v_ports)
                        for dpindex in range(len(dports)):
                            dn1, dn2 = dports[dpindex]
                            if n1:
                                if dn1:
                                    J[index * n_of_var + n1 - 1, index * n_of_var + dn1 - 1] = \
                                        J[index * n_of_var + n1 - 1, index * n_of_var + dn1 - 1] + \
                                            elem.g(opindex, v_ports, dpindex)
                                if dn2:
                                    J[index * n_of_var + n1 - 1, index * n_of_var + dn2 - 1] =\
                                        J[index * n_of_var + n1 - 1, index * n_of_var + dn2 - 1] - 1.0 * \
                                            elem.g(opindex, v_ports, dpindex)
                            if n2:
                                if dn1:
                                    J[index * n_of_var + n2 - 1, index * n_of_var + dn1 - 1] = \
                                        J[index * n_of_var + n2 - 1, index * n_of_var + dn1 - 1] - 1.0 * \
                                            elem.g(opindex, v_ports, dpindex)
                                if dn2:
                                    J[index * n_of_var + n2 - 1, index * n_of_var + dn2 - 1] =\
                                        J[index * n_of_var + n2 - 1, index * n_of_var + dn2 - 1] + \
                                            elem.g(opindex, v_ports, dpindex)

        J = J + CMAT
        residuo = CMAT*x + T + Tf + Tt
        if sparse:
            J = scipy.sparse.csc_matrix(J)
            lu = scipy.sparse.linalg.splu(J)
            dx = lu.solve(-residuo)
        else:
            dx = np.linalg.solve(J, -residuo)
        # td
        for index in range(points):
            td[index, 0] = dc_analysis.get_td(dx[index*n_of_var:(index + 1)*n_of_var,
                                                 0].reshape(-1, 1),
                                              locked_nodes, n=-1)
        x = x + min(abs(td))[0] * dx
        # convergence check
        converged = _convergence_check(dx, x, nv_indices, ni_indices,
                                      vector_norm)
        if converged:
            break
        tick.step()

        if options.bfpss_max_nr_iter and iteration == options.bfpss_max_nr_iter:
            printing.print_general_error("Hit BFPSS_MAX_NR_ITER (" +
                                         str(options.bfpss_max_nr_iter) +
                                         "), iteration halted.")
            converged = False
            break
        else:
            iteration = iteration + 1

    tick.hide(verbose > 2)
    if converged:
        printing.print_info_line(("done.", 3), verbose)
        t = np.arange(points) * step
        t = t.reshape((1, points))
        x = x.reshape((points, n_of_var))
        sol = results.pss_solution(circ=circ, method="brute-force",
                                   period=period, outfile=outfile)
        sol.set_results(t, x.T)
    else:
        print("failed.")
        sol = None
    return sol


def _convergence_check(dx, x, nv_indices, ni_indices, vector_norm):
    """Perform a convergence check using the specified vector norm"""
    # sometimes something diverges... look out
    if vector_norm(dx) is np.nan:
        raise OverflowError
    dxc = np.array(dx)
    xc = np.array(x)
    ret = (vector_norm(dxc[nv_indices]) < options.ver * \
           vector_norm(xc[nv_indices]) + options.vea) and \
          (not len(ni_indices) or \
           vector_norm(dxc[ni_indices]) < options.ier * \
           vector_norm(xc[ni_indices]) + options.iea)
    return ret


def _build_CMAT(mna, D, step, points, tick, n_of_var=None, sparse=False, verbose=3):
    if n_of_var is None:
        n_of_var = mna.shape[0]
    printing.print_info_line(("Building CMAT (%dx%d)... " %
                             (n_of_var*points, n_of_var*points), 5), verbose, print_nl=False)
    tick.reset()
    tick.display(verbose > 2)
    C1, C0 = implicit_euler.get_df_coeff(step)
    I = np.eye(n_of_var)
    M = mna + C1*D
    N = C0*D
    if sparse:
        CMAT = scipy.sparse.lil_matrix((n_of_var*points, n_of_var*points))
    else:
        CMAT = np.zeros((n_of_var*points, n_of_var*points))
    for li in range(points):  # li = line index
        for ci in range(points):
            if li == 0:
                if ci == 0:
                    temp = I
                elif ci == points - 1:
                    temp = -I
                else:
                    continue  # temp = Z
            else:
                if ci == li:
                    temp = M
                elif ci == li - 1:
                    temp = N
                else:
                    continue  # temp = Z
            CMAT = utilities.set_submatrix(row=li*n_of_var, col=ci*n_of_var,
                                           dest_matrix=CMAT, source_matrix=temp)
        tick.step()
    tick.hide(verbose > 2)
    if sparse:
        CMAT = scipy.sparse.coo_matrix(CMAT)
    printing.print_info_line(("done.", 5), verbose)

    return CMAT


def _build_x(mna, step, points, tick, x0=None, n_of_var=None, verbose=3):
    if n_of_var is None:
        n_of_var = mna.shape[0]
    printing.print_info_line(("Building x...", 5), verbose, print_nl=False)
    tick.reset()
    tick.display(verbose > 2)
    x = np.zeros((points * n_of_var, 1))
    if x0 is not None:
        if isinstance(x0, results.op_solution):
            x0 = x0.asarray()
        if x0.shape[0] != n_of_var:
            print("Warning x0 has the wrong dimensions. Using all 0s.")
        else:
            for index in range(points):
                x = utilities.set_submatrix(row=index*n_of_var, col=0,
                                            dest_matrix=x, source_matrix=x0)
                tick.step()

    tick.hide(verbose > 2)
    printing.print_info_line(("done.", 5), verbose)
    return x


def _build_Tf(sTf, points, tick, n_of_var, verbose=3):
    printing.print_info_line(("Building Tf...", 5), verbose, print_nl=False)
    tick.reset()
    tick.display(verbose > 2)
    Tf = np.zeros((points * n_of_var, 1))

    for index in range(1, points):
        Tf = utilities.set_submatrix(row=index*n_of_var, col=0, dest_matrix=Tf,
                                     source_matrix=sTf)
        tick.step()
    tick.hide(verbose > 2)
    printing.print_info_line(("done.", 5), verbose)
    return Tf


def _build_Tt(circ, points, step, tick, n_of_var, verbose=3):
    nv = circ.get_nodes_number()
    printing.print_info_line(("Building Tt...", 5), verbose, print_nl=False)
    tick.reset()
    tick.display(verbose > 2)
    Tt = np.zeros((points * n_of_var, 1))
    for index in range(1, points):
        v_eq = 0
        time = index * step
        for elem in circ:
            if (isinstance(elem, devices.VSource) or isinstance(elem, devices.ISource)) and elem.is_timedependent:
                if isinstance(elem, devices.VSource):
                    Tt[index * n_of_var + nv - 1 + v_eq, 0] = - \
                        1.0 * elem.V(time)
                elif isinstance(elem, devices.ISource):
                    if elem.n1:
                        Tt[index * n_of_var + elem.n1 - 1, 0] = \
                            Tt[index * n_of_var + elem.n1 - 1, 0] + \
                            elem.I(time)
                    if elem.n2:
                        Tt[index * n_of_var + elem.n2 - 1, 0] = \
                            Tt[index * n_of_var + elem.n2 - 1, 0] - \
                                elem.I(time)
            if circuit.is_elem_voltage_defined(elem):
                v_eq = v_eq + 1
            # print Tt[index*n_of_var:(index+1)*n_of_var]
        tick.step()
    tick.hide(verbose > 2)
    printing.print_info_line(("done.", 5), verbose)

    return Tt
