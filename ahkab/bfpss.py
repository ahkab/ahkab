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

__version__ = "0.091"

import sys
import numpy
import transient
import implicit_euler
import dc_analysis
import ticker
import options
import circuit
import printing
import utilities
import results
import devices


def bfpss(circ, period, step=None, points=None, autonomous=False, x0=None,
          mna=None, Tf=None, D=None, outfile='stdout', vector_norm=lambda v: max(abs(v)), verbose=3):
    """Performs a PSS analysis.

    Time step is constant, IE will be used as DF

    Parameters:
    circ is the circuit description class
    period is the period of the solution
    mna, D, Tf are not compulsory they will be computed if they're set to None
    step is the time step between consecutive points
    points is the number of points to be used
    step and points are mutually exclusive options:
    - if step is specified, the number of points will be automatically determined
    - if points is set, the step will be automatically determined
    - if none of them is set, options.shooting_default_points will be used as points
    autonomous has to be False, autonomous circuits are not supported
    x0 is the initial guess to be used. Needs work.
    outfile is the output filename. Defaults to stdout.
    verbose is set to zero (print errors only) if datafilename == 'stdout'.

    Returns: nothing
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
        sys.exit(0)

    (points, step) = check_step_and_points(step, points, period)

    n_of_var = mna.shape[0]
    locked_nodes = circ.get_locked_nodes()
    tick = ticker.ticker(increments_for_step=1)

    CMAT = build_CMAT(mna, D, step, points, tick, n_of_var=n_of_var,
                      verbose=verbose)

    x = build_x(mna, step, points, tick, x0=x0, n_of_var=n_of_var,
                verbose=verbose)

    Tf = build_Tf(Tf, points, tick, n_of_var=n_of_var, verbose=verbose)

    # time variable component: Tt this is always the same in each iter. So we build it once for all
    # this holds all time-dependent sources (both V/I).
    Tt = build_Tt(circ, points, step, tick, n_of_var=n_of_var, verbose=verbose)

    # Indices to differentiate between currents and voltages in the
    # convergence check
    nv_indices = []
    ni_indices = []
    nv_1 = len(circ.nodes_dict) - 1
    ni = n_of_var - nv_1
    for i in range(points):
        nv_indices += (i * mna.shape[0] * numpy.ones(nv_1) + \
                      numpy.arange(nv_1)).tolist()
        ni_indices += (i * mna.shape[0] * numpy.ones(ni) + \
                      numpy.arange(nv_1, n_of_var)).tolist()

    converged = False

    printing.print_info_line(("Solving... ", 3), verbose, print_nl=False)
    tick.reset()
    tick.display(verbose > 2)
    J = numpy.mat(numpy.zeros(CMAT.shape))
    T = numpy.mat(numpy.zeros((CMAT.shape[0], 1)))
    # td is a numpy matrix that will hold the damping factors
    td = numpy.mat(numpy.zeros((points, 1)))
    iteration = 0  # newton iteration counter

    while True:
        if iteration:  # the first time are already all zeros
            J[:, :] = 0
            T[:, 0] = 0
            td[:, 0] = 0
        for index in xrange(1, points):
            for elem in circ:
                # build all dT(xn)/dxn (stored in J) and T(x)
                if elem.is_nonlinear:
                    oports = elem.get_output_ports()
                    for opindex in range(len(oports)):
                        dports = elem.get_drive_ports(opindex)
                        v_ports = []
                        for dpindex in range(len(dports)):
                            dn1, dn2 = dports[dpindex]
                            v = 0  # build v: remember we trashed the 0 row and 0 col of mna -> -1
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
        residuo = CMAT * x + T + Tf + Tt
        dx = -1 * (numpy.linalg.inv(J) * residuo)
        # td
        for index in xrange(points):
            td[index, 0] = dc_analysis.get_td(
                dx[index * n_of_var:(index + 1) * n_of_var, 0], locked_nodes, n=-1)
        x = x + min(abs(td))[0, 0] * dx
        # convergence check
        converged = convergence_check(
            dx, x, nv_indices, ni_indices, vector_norm)
        if converged:
            break
        tick.step(verbose > 2)

        if options.shooting_max_nr_iter and iteration == options.shooting_max_nr_iter:
            printing.print_general_error(
                "Hitted SHOOTING_MAX_NR_ITER (" + str(options.shooting_max_nr_iter) + "), iteration halted.")
            converged = False
            break
        else:
            iteration = iteration + 1

    tick.hide(verbose > 2)
    if converged:
        printing.print_info_line(("done.", 3), verbose)
        t = numpy.mat(numpy.arange(points) * step)
        t = t.reshape((1, points))
        x = x.reshape((points, n_of_var))
        sol = results.pss_solution(
            circ=circ, method="brute-force", period=period, outfile=outfile, t_array=t, x_array=x.T)
    else:
        print "failed."
        sol = None
    return sol


def convergence_check(dx, x, nv_indices, ni_indices, vector_norm):
    if vector_norm(dx) is numpy.nan:  # sometimes something diverges... look out
        raise OverflowError
    dxc = numpy.array(dx)
    xc = numpy.array(x)
    ret = (vector_norm(dxc[nv_indices]) < options.ver * vector_norm(xc[nv_indices]) + options.vea) and \
          (not len(ni_indices) or \
           vector_norm(dxc[ni_indices]) < options.ier * vector_norm(xc[ni_indices]) + options.iea)
    return ret


def set_submatrix(row, col, dest_matrix, source_matrix):
    """Copies the source_matrix in dest_matrix,
    the position of the upper left corner of source matrix is (row, col) within dest_matrix

    Returns dest_matrix
    """
    for li in xrange(source_matrix.shape[0]):
        for ci in xrange(source_matrix.shape[1]):
            if source_matrix[li, ci] != 0:
                dest_matrix[row + li, col + ci] = source_matrix[li, ci]
    return dest_matrix


def get_e(index, length):
    """Builds a e(j=index) col vector
    e(j) is defined as:
        e(i, 0) = 1 if i=j
        e(i, 0) = 0 otherwise

    Returns: e(index)
    """
    e = numpy.mat(numpy.zeros((length, 1)))
    e[index, 0] = 1
    return e


def check_step_and_points(step, points, period):
    """Sets consistently the step size and the number of points, according to the given period
    Returns: (points, step)
    """
    if step is None and points is None:
        print "Warning: shooting had no step nor n. of points setted. Using", options.shooting_default_points, "points."
        points = options.shooting_default_points
    elif step is not None and points is not None:
        print "Warning: shooting had both step and n. of points setted. Using", step, "step. (NA)"
        points = None

    if points:
        step = (1.0 * period) / (points - 1)
    else:
        points = (1.0 * period) / step
        if points % 1 != 0:
            step = step + (step * (points % 1)) / int(points)
            points = int((1.0 * period) / step)
            printing.print_warning("adapted step is %g" % (step,))
        else:
            points = int(points)
        points = points + \
            1  # 0 - N where xN is in reality the first point of the second period!!

    return (points, step)


def build_CMAT(mna, D, step, points, tick, n_of_var=None, verbose=3):
    if n_of_var is None:
        n_of_var = mna.shape[0]
    printing.print_info_line(("Building CMAT (%dx%d)... " %
                             (n_of_var * points, n_of_var * points), 5), verbose, print_nl=False)
    tick.reset()
    tick.display(verbose > 2)
    (C1, C0) = implicit_euler.get_df_coeff(step)
    I = numpy.mat(numpy.eye(n_of_var))
    M = mna + C1 * D
    N = C0 * D
    # Z = numpy.mat(numpy.zeros((n_of_var, n_of_var)))
    CMAT = numpy.mat(numpy.zeros((n_of_var * points, n_of_var * points)))
    for li in xrange(points):  # li = line index
        for ci in xrange(points):
            if li == 0:
                if ci == 0:
                    temp = 1.0 * I
                elif ci == points - 1:
                    temp = -1.0 * I
                else:
                    continue  # temp = Z
            else:
                if ci == li:
                    temp = M
                elif ci == li - 1:
                    temp = N
                else:
                    continue  # temp = Z
            CMAT = set_submatrix(
                row=li * n_of_var, col=ci * n_of_var, dest_matrix=CMAT, source_matrix=temp)
        tick.step(verbose > 2)
    tick.hide(verbose > 2)
    printing.print_info_line(("done.", 5), verbose)

    return CMAT


def build_x(mna, step, points, tick, x0=None, n_of_var=None, verbose=3):
    if n_of_var is None:
        n_of_var = mna.shape[0]
    printing.print_info_line(("Building x...", 5), verbose, print_nl=False)
    tick.reset()
    tick.display(verbose > 2)
    x = numpy.mat(numpy.zeros((points * n_of_var, 1)))
    if x0 is not None:
        if isinstance(x0, results.op_solution):
            x0 = x0.asmatrix()
        if x0.shape[0] != n_of_var:
            print "Warning x0 has the wrong dimensions. Using all 0s."
        else:
            for index in xrange(points):
                x = set_submatrix(
                    row=index * n_of_var, col=0, dest_matrix=x, source_matrix=x0)
                tick.step(verbose > 2)

    tick.hide(verbose > 2)
    printing.print_info_line(("done.", 5), verbose)

    return x


def build_Tf(sTf, points, tick, n_of_var, verbose=3):
    printing.print_info_line(("Building Tf...", 5), verbose, print_nl=False)
    tick.reset()
    tick.display(verbose > 2)
    Tf = numpy.mat(numpy.zeros((points * n_of_var, 1)))

    for index in xrange(1, points):
        Tf = set_submatrix(
            row=index * n_of_var, col=0, dest_matrix=Tf, source_matrix=sTf)
        tick.step(verbose > 2)

    tick.hide(verbose > 2)
    printing.print_info_line(("done.", 5), verbose)

    return Tf


def build_Tt(circ, points, step, tick, n_of_var, verbose=3):
    nv = len(circ.nodes_dict)
    printing.print_info_line(("Building Tt...", 5), verbose, print_nl=False)
    tick.reset()
    tick.display(verbose > 2)
    Tt = numpy.zeros((points * n_of_var, 1))
    for index in xrange(1, points):
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
        tick.step(verbose > 2)
    tick.hide(verbose > 2)
    printing.print_info_line(("done.", 5), verbose)

    return Tt
