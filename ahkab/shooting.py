# -*- coding: iso-8859-1 -*-
# shooting.py
# Shooting analysis module
# Copyright 2009 Giuseppe Venturini

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

"""Periodic steady state analysis based on the shooting method."""

from __future__ import (unicode_literals, absolute_import,
                        division, print_function)

import numpy as np
import numpy.linalg

from . import transient
from . import implicit_euler
from . import dc_analysis
from . import ticker
from . import options
from . import circuit
from . import printing
from . import utilities
from . import results
from . import devices


def shooting_analysis(circ, period, step=None, x0=None, points=None, autonomous=False,
             matrices=None, outfile='stdout', vector_norm=lambda v: max(abs(v)), verbose=3):
    """Performs a periodic steady state analysis based on the algorithm described in:

        Brambilla, A.; D'Amore, D., "Method for steady-state simulation of
        strongly nonlinear circuits in the time domain," *Circuits and
        Systems I: Fundamental Theory and Applications, IEEE Transactions on*,
        vol.48, no.7, pp.885-889, Jul 2001.

        http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=933329&isnumber=20194

    The results have been computed again by me, the formulation is not exactly the
    same, but the idea behind the shooting algorithm is.

    This method allows us to have a period with many points without having to
    invert a huge matrix (and being limited to the maximum matrix size).

    A transient analysis is performed to initialize the solver.

    We compute the change in the last point, calculating several matrices in
    the process.
    From that, with the same matrices we calculate the changes in all points,
    starting from 0 (which is the same as the last one), then 1, ...

    Key points:

    - Only non-autonomous circuits are supported.
    - The time step is constant.
    - Implicit Euler is used as DF.

    **Parameters:**

    circ : Circuit instance
        The circuit description class.
    period : float
        The period of the solution.
    step : float, optional
        The time step between consecutive points.
        If not set, it will be computed from ``period`` and ``points``.
    points : int, optional
        The number of points to be used. If not set, it will be computed from
        ``period`` and ``step``.
    autonomous : bool, optional
        This parameter has to be ``False``, autonomous circuits are not
        currently supported.
    matrices : dict, optional
        A dictionary that may have as keys 'MNA', 'N' and 'D', with entries set
        to the corresponding MNA-formulation matrices, in case they have been
        already computed and the user wishes to save time by reusing them.
        Defaults to ``None`` (recompute).
    outfile : string, optional
        The output filename. Please use ``stdout`` (the default) to print to the
        standard output.
    verbose : boolean, optional
        Verbosity switch (0-6). It is set to zero (print errors only)
        if ``outfile`` == 'stdout'``, as not to corrupt the data.

    Notice that ``step`` and ``points`` are mutually exclusive options:

    - if ``step`` is specified, the number of points will be automatically determined.
    - if ``points`` is set, the step will be automatically determined.
    - if none of them is set, ``options.shooting_default_points`` will be used as points.

    **Returns:**

    sol : PSS solution object or ``None``
        The solution. If the circuit can't be solve, ``None`` is returned instead.
    """

    if outfile == "stdout":
        verbose = 0

    printing.print_info_line(("Starting periodic steady state analysis:", 3),
                             verbose)
    printing.print_info_line(("Method: shooting", 3), verbose)

    if isinstance(x0, results.op_solution):
        x0 = x0.asarray()

    # Do we have MNA T D or do we need to build them
    # MNA & T
    if (matrices is None or type(matrices) != dict or 'MNA' not in matrices or
        'Tf' not in matrices):
        # recalculate
        mna, Tf = dc_analysis.generate_mna_and_N(circ, verbose=verbose)
        mna = utilities.remove_row_and_col(mna)
        Tf = utilities.remove_row(Tf, rrow=0)
    elif not matrices['MNA'].shape[0] == matrices['Tf'].shape[0]:
        raise ValueError("MNA matrix and N vector have different number of" +
                         " rows.")
    else:
        mna, Tf = matrices['MNA'], matrices['Tf']
    # D
    if matrices is None or 'D' not in matrices or matrices['D'] is None:
        D = transient.generate_D(circ, [mna.shape[0], mna.shape[0]])
        D = utilities.remove_row_and_col(D)
    elif not mna.shape == matrices['D'].shape:
        raise ValueError("MNA matrix and D matrix have different sizes.")
    else:
        D = matrices['D']

    points, step = utilities.check_step_and_points(step, points, period,
                                                   options.shooting_default_points)

    n_of_var = mna.shape[0]
    locked_nodes = circ.get_locked_nodes()

    printing.print_info_line(("Starting TRAN analysis for algorithm init: " +
                              ("stop=%g, step=%g... " % (10*points*step, step)),
                              3), verbose, print_nl=False)
    xtran = transient.transient_analysis(circ=circ, tstart=0, tstep=step,
                                         tstop=10*points*step, method="TRAP",
                                         x0=None, mna=mna, N=Tf, D=D,
                                         use_step_control=False,
                                         outfile=outfile+".tran",
                                         return_req_dict={"points": points},
                                         verbose=0)
    if xtran is None:
        print("failed.")
        return None
    printing.print_info_line(("done.", 3), verbose)

    x = []
    for index in range(points):
        x.append(xtran[index * n_of_var:(index + 1) * n_of_var, 0])

    tick = ticker.ticker(increments_for_step=1)

    MAass_static, MBass = _build_static_MAass_and_MBass(mna, D, step)

    # This contains the time invariant part, Tf. The time variable component,
    # Tt, is always the same, since the time interval is the same this holds all
    # time-dependent sources (both V/I).

    # silly numpy sum array of size 6,1 with array of size 6
    # get array of size 6,6
    Tf = Tf.squeeze()
    Tass_static_vector = _build_Tass_static_vector(circ, Tf, points, step, tick,
                                                   n_of_var, verbose)

    converged = False
    printing.print_info_line(("Solving... ", 3), verbose, print_nl=False)
    tick.reset()
    tick.display(verbose > 2)

    iteration = 0  # newton iteration counter
    conv_counter = 0

    while True:
        dx = []
        Tass_variable_vector = []
        MAass_variable_vector = []
        for index in range(points):
            if index == 0:
                xn_minus_1 = x[points - 1]
            else:
                xn_minus_1 = x[index - 1]
            MAass_variable, Tass_variable = _get_variable_MAass_and_Tass(circ,
                    x[index], xn_minus_1, mna, D, step, n_of_var)
            MAass_variable_vector.append(MAass_variable + MAass_static)
            Tass_variable_vector.append(Tass_variable +
                                        Tass_static_vector[index])

        dxN = _compute_dxN(MAass_variable_vector, MBass, Tass_variable_vector,
                           n_of_var, points)
        td = dc_analysis.get_td(dxN.reshape(-1, 1), locked_nodes, n=-1)
        x[points - 1] = td * dxN + x[points - 1]

        for index in range(points - 1):
            if index == 0:
                dxi_minus_1 = dxN
            else:
                dxi_minus_1 = dx[index - 1]
            dx.append(_compute_dx(MAass_variable_vector[index], MBass,
                      Tass_variable_vector[index], dxi_minus_1))
            td = dc_analysis.get_td(dx[index].reshape(-1, 1), locked_nodes,
                                    n=-1)
            x[index] = td*dx[index] + x[index]
        dx.append(dxN)

        if (_vector_norm_wrapper(dx, vector_norm) < min(options.ver, options.ier) *
                _vector_norm_wrapper(x, vector_norm) + min(options.vea, options.iea)):
            # and (dc_analysis.vector_norm(residuo) <
            # options.er*dc_analysis.vector_norm(x) + options.ea):
            if conv_counter == 3:
                converged = True
                break
            else:
                conv_counter = conv_counter + 1
        elif vector_norm(dx[points - 1]) is np.nan:  # needs work fixme
            raise OverflowError
            # break
        else:
            conv_counter = 0
            tick.step()

        if options.shooting_max_nr_iter and iteration == options.shooting_max_nr_iter:
            printing.print_general_error("Hitted SHOOTING_MAX_NR_ITER (" +
                                         str(options.shooting_max_nr_iter) +
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
        xmat = x[0].reshape(-1, 1)
        for index in range(1, points):
            xmat = numpy.concatenate((xmat, x[index].reshape(-1, 1)), axis=1)
        sol = results.pss_solution(circ=circ, method="shooting", period=period,
                                   outfile=outfile)
        sol.set_results(t, xmat)
        # print_results(circ, x, fdata, points, step)
    else:
        print("failed.")
        sol = None
    return sol


def _vector_norm_wrapper(vector, norm_fun):
    _max = 0
    for elem in vector:
        new_max = norm_fun(elem)
        if _max < new_max:
            _max = new_max
    return _max


def _build_static_MAass_and_MBass(mna, D, step):
    (C1, C0) = implicit_euler.get_df_coeff(step)
    MAass = mna + D * C1
    MBass = np.dot(D, C0)
    return (MAass, MBass)


def _build_Tass_static_vector(circ, Tf, points, step, tick, n_of_var, verbose=3):
    Tass_vector = []
    nv = circ.get_nodes_number()
    printing.print_info_line(("Building Tass...", 5), verbose, print_nl=False)

    tick.reset()
    tick.display(verbose > 2)
    for index in range(0, points):
        Tt = numpy.zeros((n_of_var,))
        v_eq = 0
        time = index * step
        for elem in circ:
            if (isinstance(elem, devices.VSource) or
                isinstance(elem, devices.ISource)) and elem.is_timedependent:
                # time dependent source
                if isinstance(elem, devices.VSource):
                    Tt[nv - 1 + v_eq] = -1.0 * elem.V(time)
                elif isinstance(elem, devices.ISource):
                    if elem.n1:
                        Tt[elem.n1 - 1] = Tt[elem.n1 - 1] + elem.I(time)
                    if elem.n2:
                        Tt[elem.n2 - 1] = Tt[elem.n2 - 1] - elem.I(time)
            if circuit.is_elem_voltage_defined(elem):
                v_eq = v_eq + 1
        tick.step()
        Tass_vector.append(Tf + Tt)
    tick.hide(verbose > 2)
    printing.print_info_line(("done.", 5), verbose)

    return Tass_vector


def _get_variable_MAass_and_Tass(circ, xi, xi_minus_1, M, D, step, n_of_var):
    Tass = np.zeros((n_of_var,))
    J = np.zeros((n_of_var, n_of_var))
    (C1, C0) = implicit_euler.get_df_coeff(step)

    for elem in circ:
        # build all dT(xn)/dxn (stored in J) and T(x)
        if elem.is_nonlinear:
            output_ports = elem.get_output_ports()
            for index in range(len(output_ports)):
                n1, n2 = output_ports[index]
                ports = elem.get_drive_ports(index)
                v_ports = []
                for port in ports:
                    v = 0  # build v: remember we trashed the 0 row and 0 col of mna -> -1
                    if port[0]:
                        v = v + xi[port[0] - 1]
                    if port[1]:
                        v = v - xi[port[1] - 1]
                    v_ports.append(v)
                if n1:
                    Tass[n1 - 1] = Tass[n1 - 1] + elem.i(index, v_ports)
                if n2:
                    Tass[n2 - 1] = Tass[n2 - 1] - elem.i(index, v_ports)
                for pindex in range(len(ports)):
                    if n1:
                        if ports[pindex][0]:
                            J[n1 - 1, ports[pindex][0] - 1] = \
                                J[n1 - 1, ports[pindex][0] - 1] + \
                                elem.g(index, v_ports, pindex)
                        if ports[pindex][1]:
                            J[n1 - 1, ports[pindex][1] - 1] = \
                                J[n1 - 1, ports[pindex][1] - 1] - \
                                elem.g(index, v_ports, pindex)
                    if n2:
                        if ports[pindex][0]:
                            J[n2 - 1, ports[pindex][0] - 1] = \
                                J[n2 - 1, ports[pindex][0] - 1] - \
                                elem.g(index, v_ports, pindex)
                        if ports[pindex][1]:
                            J[n2 - 1, ports[pindex][1] - 1] = \
                                J[n2 - 1, ports[pindex][1] - 1] + \
                                elem.g(index, v_ports, pindex)

    Tass = Tass + np.dot(D*C1, xi) + np.dot(M, xi) + np.dot(D * C0, xi_minus_1)

    return J, Tass


def _compute_dxN(MAass_vector, MBass, Tass_vector, n_of_var, points):
    temp_mat1 = np.eye(n_of_var)
    for index in range(points):
        temp_mat1 = -np.linalg.solve(MAass_vector[index], np.dot(MBass, temp_mat1))
    temp_mat2 = np.zeros((n_of_var,))
    for index in range(points):
        temp_mat3 = -np.linalg.solve(MAass_vector[index], Tass_vector[index])
        for index2 in range(index + 1, points):
            temp_mat3 = -np.linalg.solve(MAass_vector[index2], np.dot(MBass, temp_mat3))
        temp_mat2 = temp_mat2 + temp_mat3

    dxN = np.dot(np.linalg.inv(np.eye(n_of_var) - temp_mat1), temp_mat2)

    return dxN


def _compute_dx(MAass, MBass, Tass, dxi_minus_1):
    dxi = -np.linalg.solve(MAass, np.dot(MBass, dxi_minus_1) + Tass)
    # dxi = -1 * np.linalg.inv(MAass) * (MBass * dxi_minus_1 + Tass)
    return dxi
