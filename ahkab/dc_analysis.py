# -*- coding: iso-8859-1 -*-
# dc_analysis.py
# DC simulation methods
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

"""
This module provides the functions needed to perform OP and DC simulations.

The principal are:

* :func:`dc_analysis` - which performs a dc sweep,
* :func:`op_analysis` - which does an operation point analysis or

Notice that internally, :func:`dc_analysis` calls :func:`op_analysis`,
since a DC sweep is nothing but a series of OP analyses..

The actual circuit solution is done by :func:`mdn_solver`, that uses a
modified version of the Newton Rhapson method.

Module reference
################

"""

from __future__ import (unicode_literals, absolute_import,
                        division, print_function)

import sys
import re
import copy

import numpy as np
import numpy.linalg
import scipy.sparse
import scipy.sparse.linalg

from . import devices
from . import diode
from . import constants
from . import ticker
from . import options
from . import circuit
from . import printing
from . import utilities
from . import dc_guess
from . import results

from .utilities import convergence_check

specs = {'op': {
    'tokens': ({
               'label': 'guess',
               'pos': None,
               'type': bool,
               'needed': False,
               'dest': 'guess',
               'default': options.dc_use_guess
               },
        {
               'label': 'ic_label',
               'pos': None,
               'type': str,
               'needed': False,
               'dest': 'x0',
               'default': None
               }
               )
},
    'dc': {'tokens': ({
                      'label': 'source',
                      'pos': 0,
                      'type': str,
                      'needed': True,
                      'dest': 'source',
                      'default': None
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
                      'label': 'step',
                      'pos': 3,
                      'type': float,
                      'needed': True,
                      'dest': 'step',
                      'default': None
                      },
                      {
                      'label': 'type',
                      'pos': None,
                      'type': str,
                      'needed': False,
                      'dest': 'sweep_type',
                      'default': options.dc_lin_step
                      }
                     )
           }
}


def dc_solve(mna, Ndc, circ, Ntran=None, Gmin=None, x0=None, time=None,
             MAXIT=None, locked_nodes=None, skip_Tt=False, verbose=3):
    """Low-level method to perform a DC solution of the circuit

    .. note::
        Typically the user calls :func:`dc_analysis.op_analysis` or
        :func:`dc_analysis.dc_analysis`, which in turn will setup all
        matrices and call this method on their behalf.

    The system we want to solve is:

    .. math::

        (mna + G_{min}) \\cdot x + N(t) + T(x, t) = 0

    Where:

    * :math:`mna` is the reduced MNA matrix with the required KVL/KCL rows
    * :math:`N` is composed by a DC part, :math:`N_{dc}`, and a dynamic
      time-dependent part :math:`N_{tran}(t)` and a time-dependent part
      :math:`T_t(t)`.
    * :math:`T(x, t)` is both time-dependent and non-linear with respect to
      the circuit solution :math:`x`, and it will be built at each iteration
      over :math:`t` and :math:`x`.

    **Parameters:**

    mna : ndarray
        The MNA matrix described above. It can be built calling
        :func:`generate_mna_and_N`. This matrix will contain the dynamic
        component due to a Differetiation Formula (DF) when this method is
        called from a transient analysis.
    Ndc : ndarray
        The DC part of :math:`N`. Also this vector may be built calling
        :func:`generate_mna_and_N`.
    circ : Circuit instance
        The circuit instance from which ``mna`` and ``N`` were built.
    Ntran : ndarray, optional
        The linear time-dependent and *dynamic* part of :math:`N`, if available.
        Notice this is typically set when a DF being applied and the method is
        being called from a transient analysis.
    Gmin : ndarray, optional
        A matrix of the same size of ``mna``, containing the minimum
        transconductances to ground. It can be built with
        :func:`build_gmin_matrix`. If not set, no Gmin matrix is used.
    x0 : ndarray or results.op_solution instance, optional
        The initial guess for the Newthon-Rhapson algorithm. If not specified,
        the all-zeros vector will be used.
    time : float scalar, optional
        The time at which any matrix evaluation done by this method will be
        performed. Do not set for DC or OP analysis, must be set for a
        transisent analysis. Notice that :math:`t=0` is not the same as DC!
    MAXIT : int, optional
        The maximum number of Newton Rhapson iterations to be performed before
        giving up. If unset, ``options.dc_max_nr_iter`` is used.
    locked_nodes : list of tuples, optional
        The nodes that need to have a well behaved, slowly varying voltage
        applied. Typically they control non-linear elements. This is generated
        by :func:`ahkab.circuit.Circuit.get_locked_nodes` and it will be
        generated for you if left unset. However, if you are doing many
        simulations of the same circuit (as it happens in a transient
        analysis), it's a good idea to generate it only once.
    skip_Tt : boolean, optional
        Do not build the :math:`T_t(t)` vector. Defaults to ``False``.
    verbose : int, optional
        The verbosity level. From 0 (silent) to 6 (debug). Defaults to 3.

    **Returns:**

    x : ndarray
        The solution, if found.
    error : ndarray
        The error associated with each solution item, if it was found.
    converged : boolean
        A flag set to True when convergence was detected.
    tot_iterations : int
        Total number of NR iterations run.
    """
    if MAXIT == None:
        MAXIT = options.dc_max_nr_iter
    if locked_nodes is None:
        locked_nodes = circ.get_locked_nodes()
    mna_size = mna.shape[0]
    nv = circ.get_nodes_number()
    tot_iterations = 0

    if Gmin is None:
        Gmin = 0

    if Ntran is None:
        Ntran = 0

    # time variable component: Tt this is always the same in each iter. So we
    # build it once for all.
    Tt = np.zeros((mna_size, 1))
    v_eq = 0
    if not skip_Tt:
        for elem in circ:
            if (isinstance(elem, devices.VSource) or isinstance(elem, devices.ISource)) and elem.is_timedependent:
                if isinstance(elem, devices.VSource):
                    Tt[nv - 1 + v_eq, 0] = -1 * elem.V(time)
                elif isinstance(elem, devices.ISource):
                    if elem.n1:
                        Tt[elem.n1 - 1, 0] = Tt[elem.n1 - 1, 0] + elem.I(time)
                    if elem.n2:
                        Tt[elem.n2 - 1, 0] = Tt[elem.n2 - 1, 0] - elem.I(time)
            if circuit.is_elem_voltage_defined(elem):
                v_eq = v_eq + 1
    # update N to include the time variable sources
    Ndc = Ndc + Tt

    # initial guess, if specified, otherwise it's zero
    if x0 is not None:
        if isinstance(x0, results.op_solution):
            x = x0.asarray()
        else:
            x = x0
    else:
        x = np.zeros((mna_size, 1))
                      # has n-1 rows because of discard of ^^^

    converged = False
    standard_solving, gmin_stepping, source_stepping = get_solve_methods()
    standard_solving, gmin_stepping, source_stepping = set_next_solve_method(
        standard_solving, gmin_stepping,
        source_stepping, verbose)

    convergence_by_node = None
    printing.print_info_line(("Solving... ", 3), verbose, print_nl=False)

    while(not converged):
        if standard_solving["enabled"]:
            mna_to_pass = mna + Gmin
            N_to_pass = Ndc + Ntran * (Ntran is not None)
        elif gmin_stepping["enabled"]:
            # print "gmin index:", str(gmin_stepping["index"])+", gmin:", str(
            # 10**(gmin_stepping["factors"][gmin_stepping["index"]]))
            printing.print_info_line(
                ("Setting Gmin to: " + str(10 ** gmin_stepping["factors"][gmin_stepping["index"]]), 6), verbose)
            mna_to_pass = build_gmin_matrix(
                circ, 10 ** (gmin_stepping["factors"][gmin_stepping["index"]]), mna_size, verbose) + mna
            N_to_pass = Ndc + Ntran * (Ntran is not None)
        elif source_stepping["enabled"]:
            printing.print_info_line(
                ("Setting sources to " + str(source_stepping["factors"][source_stepping["index"]] * 100) + "% of their actual value", 6), verbose)
            mna_to_pass = mna + Gmin
            N_to_pass = source_stepping["factors"][source_stepping["index"]]*Ndc + Ntran*(Ntran is not None)
        try:
            (x, error, converged, n_iter, convergence_by_node) = mdn_solver(x, mna_to_pass, circ, T=N_to_pass,
                                                                            nv=nv, print_steps=(verbose > 0), locked_nodes=locked_nodes, time=time, MAXIT=MAXIT, debug=(verbose == 6))
            tot_iterations += n_iter
        except np.linalg.linalg.LinAlgError:
            n_iter = 0
            converged = False
            print("failed.")
            printing.print_general_error("J Matrix is singular")
        except OverflowError:
            n_iter = 0
            converged = False
            print("failed.")
            printing.print_general_error("Overflow")

        if not converged:
            if verbose == 6 and convergence_by_node is not None:
                for ivalue in range(len(convergence_by_node)):
                    if not convergence_by_node[ivalue] and ivalue < nv - 1:
                        print("Convergence problem node %s" % (circ.int_node_to_ext(ivalue),))
                    elif not convergence_by_node[ivalue] and ivalue >= nv - 1:
                        e = circ.find_vde(ivalue)
                        print("Convergence problem current in %s" % e.part_id)
            if n_iter == MAXIT - 1:
                printing.print_general_error(
                    "Error: MAXIT exceeded (" + str(MAXIT) + ")")
            if more_solve_methods_available(standard_solving, gmin_stepping, source_stepping):
                standard_solving, gmin_stepping, source_stepping = set_next_solve_method(
                    standard_solving, gmin_stepping, source_stepping, verbose)
            else:
                # print "Giving up."
                x = None
                error = None
                break
        else:
            printing.print_info_line(
                ("[%d iterations]" % (n_iter,), 6), verbose)
            if (source_stepping["enabled"] and source_stepping["index"] != 9):
                converged = False
                source_stepping["index"] = source_stepping["index"] + 1
            elif (gmin_stepping["enabled"] and gmin_stepping["index"] != 9):
                gmin_stepping["index"] = gmin_stepping["index"] + 1
                converged = False
            else:
                printing.print_info_line((" done.", 3), verbose)
    return (x, error, converged, tot_iterations)


def build_gmin_matrix(circ, gmin, mna_size, verbose):
    """Build a Gmin matrix

    **Parameters:**

    circ : circuit instance
        The circuit for which the matrix is built.
    gmin : scalar float
        The value of the minimum conductance to ground to be used.
    mna_size : int
        The size of the MNA matrix associated with the GMIN matrix being built.
    verbose : int
        The verbosity level, from 0 (silent) to 6 (debug).

    **Returns:**

    Gmin : ndarray of size (mna_size, mna_size)
        The Gmin matrix itself.

    """
    printing.print_info_line(("Building Gmin matrix...", 5), verbose)
    Gmin_matrix = np.zeros((mna_size, mna_size))
    for index in range(circ.get_nodes_number() - 1):
        Gmin_matrix[index, index] = gmin
        # the three missing terms of the stample matrix go on [index,0] [0,0] [0, index] but since
        # we discarded the 0 row and 0 column, we simply don't need to add them
        # the last lines are the KVL lines, introduced by voltage sources.
        # Don't add gmin there.
    return Gmin_matrix


def set_next_solve_method(standard_solving, gmin_stepping, source_stepping, verbose=3):
    """Select the next solving method.

    We have the standard solving method and two homotopies available. The
    homotopies are :math:`G_{min}` stepping and source stepping.

    They will be selected and enabled when failures occur according to the
    options values:

    * ``options.use_standard_solve_method``,
    * ``options.use_gmin_stepping``,
    * ``options.use_source_stepping``.

    The methods will be used in the order above.

    The inputs to this method are three dictionaries that keep track of which
    method is currently enabled and which ones has failed in the past.

    **Parameters:**

    standard_solving, gmin_stepping, source_stepping : dict
        The dictionaries contain the options and the status of the methods, they
        should be the values provided by :func:`get_solve_methods`.
    verbose : int, optional
        The verbosity level, from 0 (silent) to 6 (debug).

    **Returns:**

    standard_solving, gmin_stepping, source_stepping : dict
        The updated dictionaries.
    """
    if standard_solving["enabled"]:
        printing.print_info_line(("failed.", 1), verbose)
        standard_solving["enabled"] = False
        standard_solving["failed"] = True
    elif gmin_stepping["enabled"]:
        printing.print_info_line(("failed.", 1), verbose)
        gmin_stepping["enabled"] = False
        gmin_stepping["failed"] = True
    elif source_stepping["enabled"]:
        printing.print_info_line(("failed.", 1), verbose)
        source_stepping["enabled"] = False
        source_stepping["failed"] = True
    if not standard_solving["failed"] and options.use_standard_solve_method:
        standard_solving["enabled"] = True
    elif not gmin_stepping["failed"] and options.use_gmin_stepping:
        gmin_stepping["enabled"] = True
        printing.print_info_line(
            ("Enabling gmin stepping convergence aid.", 3), verbose)
    elif not source_stepping["failed"] and options.use_source_stepping:
        source_stepping["enabled"] = True
        printing.print_info_line(
            ("Enabling source stepping convergence aid.", 3), verbose)

    return standard_solving, gmin_stepping, source_stepping


def more_solve_methods_available(standard_solving, gmin_stepping, source_stepping):
    """Are there more solving methods available?

    **Parameters:**

    standard_solving, gmin_stepping, source_stepping : dict
        The dictionaries contain the options and the status of the methods.

    **Returns:**

    rsp : boolean
        The answer.
    """

    if (standard_solving["failed"] or not options.use_standard_solve_method) and \
       (gmin_stepping["failed"] or not options.use_gmin_stepping) and \
       (source_stepping["failed"] or not options.use_source_stepping):
        return False
    else:
        return True


def get_solve_methods():
    """Get all the available solving methods

    We have the standard solving method and two homotopies available. The
    homotopies are :math:`G_{min}` stepping and source stepping.

    Solving methods may be enabled and disabled through the
    options values:

    * ``options.use_standard_solve_method``,
    * ``options.use_gmin_stepping``,
    * ``options.use_source_stepping``.

    **Returns:**

    standard_solving, gmin_stepping, source_stepping : dict
        The dictionaries contain the options and the status of the methods.
    """
    standard_solving = {"enabled": False, "failed": False}
    g_indices = list(range(int(numpy.log(options.gmin)), 0))
    g_indices.reverse()
    gmin_stepping = {"enabled": False, "failed":
                     False, "factors": g_indices, "index": 0}
    source_stepping = {"enabled": False, "failed": False, "factors": (
        0.001, .005, .01, .03, .1, .3, .5, .7, .8, .9), "index": 0}
    return standard_solving, gmin_stepping, source_stepping


def dc_analysis(circ, start, stop, step, source, sweep_type='LINEAR', guess=True, x0=None, outfile="stdout", verbose=3):
    """Performs a sweep of the value of V or I of a independent source from start
    value to stop value using the provided step.

    For every circuit generated, computes the OP.  This function relays on
    :func:`dc_analysis.op_analysis` to actually solve each circuit.

    **Parameters:**

    circ : Circuit instance
        The circuit instance to be simulated.
    start : float
        Start value of the sweep source
    stop : float
        Stop value of the sweep source
    step : float
        The step size in the sweep
    source : string
        The part ID of the source to be swept, eg. ``'V1'``.
    sweep_type : string, optional
        Either options.dc_lin_step (default) or options.dc_log_step
    guess : boolean, optional
        op_analysis will guess to start the first NR iteration for the first point,
        the previsious OP is used from then on. Defaults to ``True``.
    outfile : string, optional
        Filename of the output file. If set to ``'stdout'`` (default), prints to
        screen.
    verbose : int
        The verbosity level, from 0 (silent) to 6 (debug).

    **Returns:**

    rstdc : results.dc_solution instance or None
        A ``results.dc_solution`` instance is returned, if a solution was found
        for at least one sweep value.  or ``None``, if an error occurred (eg
        invalid start/stop/step values) or there was no solution for any
        sweep value.
    """
    if outfile == 'stdout':
        verbose = 0
    printing.print_info_line(("Starting DC analysis:", 2), verbose)
    elem_type, elem_descr = source[0].lower(), source.lower() # eg. 'v', 'v34'
    sweep_label = elem_type[0].upper() + elem_descr[1:]

    if sweep_type == options.dc_log_step and stop - start < 0:
        printing.print_general_error(
            "DC analysis has log sweeping and negative stepping.")
        sys.exit(1)
    if (stop - start) * step < 0:
        raise ValueError("Unbonded stepping in DC analysis.")

    points = (stop - start) / step + 1
    sweep_type = sweep_type.upper()[:3]

    if sweep_type == options.dc_log_step:
        dc_iter = utilities.log_axis_iterator(start, stop, points=points)
    elif sweep_type == options.dc_lin_step:
        dc_iter = utilities.lin_axis_iterator(start, stop, points=points)
    else:
        printing.print_general_error("Unknown sweep type: %s" % (sweep_type,))
        sys.exit(1)

    if elem_type != 'v' and elem_type != 'i':
        printing.print_general_error(
            "Sweeping is possible only with voltage and current sources. (" + str(elem_type) + ")")
        sys.exit(1)

    source_elem = None
    for index in range(len(circ)):
        if circ[index].part_id.lower() == elem_descr:
            if elem_type == 'v':
                if isinstance(circ[index], devices.VSource):
                    source_elem = circ[index]
                    break
            if elem_type == 'i':
                if isinstance(circ[index], devices.ISource):
                    source_elem = circ[index]
                    break
    if not source_elem:
        raise ValueError(".DC: source %s was not found." % source)

    if isinstance(source_elem, devices.VSource):
        initial_value = source_elem.dc_value
    else:
        initial_value = source_elem.dc_value

    # If the initial value is set to None, op_analysis will attempt a smart guess (if guess),
    # Then for each iteration, the last result is used as x0, since op_analysis will not
    # attempt to guess the op if x0 is not None.
    x = x0

    sol = results.dc_solution(
        circ, start, stop, sweepvar=sweep_label, stype=sweep_type, outfile=outfile)

    printing.print_info_line(("Solving... ", 3), verbose, print_nl=False)
    tick = ticker.ticker(1)
    tick.display(verbose > 2)

    # sweep setup

    # tarocca il generatore di tensione, avvia DC silenziosa, ritarocca etc
    index = 0
    for sweep_value in dc_iter:
        index = index + 1
        if isinstance(source_elem, devices.VSource):
            source_elem.dc_value = sweep_value
        else:
            source_elem.dc_value = sweep_value
        # silently calculate the op
        x = op_analysis(circ, x0=x, guess=guess, verbose=0)
        if x is None:
            tick.hide(verbose > 2)
            if not options.dc_sweep_skip_allowed:
                print("Could't solve the circuit for sweep value:", start + index * step)
                solved = False
                break
            else:
                print("Skipping sweep value:", start + index * step)
                continue
        solved = True
        sol.add_op(sweep_value, x)

        tick.step()

    tick.hide(verbose > 2)
    if solved:
        printing.print_info_line(("done", 3), verbose)

    # clean up
    if isinstance(source_elem, devices.VSource):
        source_elem.dc_value = initial_value
    else:
        source_elem.dc_value = initial_value

    return sol if solved else None


def op_analysis(circ, x0=None, guess=True, outfile=None, verbose=3):
    """Runs an Operating Point (OP) analysis

    **Parameters:**

    circ : Circuit instance
        The circuit instance on which the simulation is run
    x0 : op_solution instance or ndarray, optional
        The initial guess to be used to start the NR :func:`mdn_solver`.
    guess : boolean, optional
        If set to ``True`` (default) and ``x0`` is ``None``, it will generate a
        'smart' guess to use as ``x0``.
    verbose : int
        The verbosity level from 0 (silent) to 6 (debug).

    **Returns:**

    A ``result.op_solution`` instance, if successful, ``None`` otherwise.
    """
    if outfile == 'stdout':
        verbose = 0  # silent mode, print out results only.
    if not options.dc_use_guess:
        guess = False

    (mna, N) = generate_mna_and_N(circ, verbose=verbose)

    printing.print_info_line(
        ("MNA matrix and constant term (complete):", 4), verbose)
    printing.print_info_line((mna, 4), verbose)
    printing.print_info_line((N, 4), verbose)

    # lets trash the unneeded col & row
    printing.print_info_line(
        ("Removing unneeded row and column...", 4), verbose)
    mna = utilities.remove_row_and_col(mna)
    N = utilities.remove_row(N, rrow=0)

    printing.print_info_line(("Starting op analysis:", 2), verbose)

    if x0 is None and guess:
        x0 = dc_guess.get_dc_guess(circ, verbose=verbose)
    # if x0 is not None, use that

    printing.print_info_line(("Solving with Gmin:", 4), verbose)
    Gmin_matrix = build_gmin_matrix(
        circ, options.gmin, mna.shape[0], verbose - 2)
    (x1, error1, solved1, n_iter1) = dc_solve(mna, N,
                                              circ, Gmin=Gmin_matrix, x0=x0, verbose=verbose)

    # We'll check the results now. Recalculate them without Gmin (using previsious solution as initial guess)
    # and check that differences on nodes and current do not exceed the
    # tolerances.
    if solved1:
        op1 = results.op_solution(
            x1, error1, circ, outfile=outfile, iterations=n_iter1)
        printing.print_info_line(("Solving without Gmin:", 4), verbose)
        (x2, error2, solved2, n_iter2) = dc_solve(
            mna, N, circ, Gmin=None, x0=x1, verbose=verbose)
    else:
        solved2 = False

    if solved1 and not solved2:
        printing.print_general_error("Can't solve without Gmin.")
        if verbose:
            print("Displaying latest valid results.")
            op1.write_to_file(filename='stdout')
        opsolution = op1
    elif solved1 and solved2:
        op2 = results.op_solution(
            x2, error2, circ, outfile=outfile, iterations=n_iter1 + n_iter2)
        op2.gmin = 0
        badvars = results.op_solution.gmin_check(op2, op1)
        printing.print_result_check(badvars, verbose=verbose)
        check_ok = not (len(badvars) > 0)
        if not check_ok and verbose:
            print("Solution with Gmin:")
            op1.write_to_file(filename='stdout')
            print("Solution without Gmin:")
            op2.write_to_file(filename='stdout')
        opsolution = op2
    else:  # not solved1
        printing.print_general_error("Couldn't solve the circuit. Giving up.")
        opsolution = None

    if opsolution and outfile != 'stdout' and outfile is not None:
        opsolution.write_to_file()
    if opsolution and (verbose > 2 or outfile == 'stdout') and options.cli:
        opsolution.write_to_file(filename='stdout')

    return opsolution


def mdn_solver(x, mna, circ, T, MAXIT, nv, locked_nodes, time=None,
               print_steps=False, vector_norm=lambda v: max(abs(v)),
               debug=True):
    """
    Solves a problem like F(x) = 0 using the Newton Algorithm with a variable
    damping.

    Where:

    .. math::

        F(x) = mna*x + T + T(x)

    * :math:`mna` is the Modified Network Analysis matrix of the circuit
    * :math:`T(x)` is the contribute of nonlinear elements to KCL
    * :math:`T` contains the contributions of the independent sources, time
    * invariant and linear

    If :math:`x(0)` is the initial guess, every :math:`x(n+1)` is given by:

    .. math::
        x(n+1) = x(n) + td \\cdot dx

    Where :math:`td` is a damping coefficient to avoid overflow in non-linear
    components and excessive oscillation in the very first iteration. Afterwards
    :math:`td=1` To calculate :math:`td`, an array of locked nodes is needed.

    The convergence check is done this way:

    **Parameters:**

    x : ndarray
        The initial guess. If set to ``None``, it will be initialized to all
        zeros. Specifying a initial guess may improve the convergence time of
        the algorithm and determine which solution (if any) is found if there
        are more than one.
    mna : ndarray
        The Modified Network Analysis matrix of the circuit, reduced, see above.
    circ : circuit instance
        The circuit instance.
    T : ndarray,
        The :math:`T` vector described above.
    MAXIT : int
        The maximum iterations that the method may perform.
    nv : int
        Number of nodes in the circuit (counting the ref, 0)
    locked_nodes : list of tuples
        A list of ports driving non-linear elements, generated by
        ``circ.get_locked_nodes()``
    time : float or None, optional
        The value of time to be passed to non_linear _and_ time variant
        elements.
    print_steps : boolean, optional
        Show a progress indicator, very verbose. Defaults to ``False``.
    vector_norm : function, optional
        An R^N -> R^1 function returning the norm of a vector, for convergence
        checking. Defaults to the maximum norm, ie :math:`f(x) = max(|x|)`,
    debug : int, optional
        Debug flag that will result in an array being returned containing
        node-by-node convergence information.

    **Returns:**

    sol : ndarray
        The solution.
    err : ndarray
        The remaining error.
    converged : boolean
        A boolean that is set to ``True`` whenever the method exits because of a
        successful convergence check. ``False`` whenever convergence problems
        where found.
    N : int
        The number of NR iterations performed.
    convergence_by_node : list
        If ``debug`` was set to ``True``, this list has the same size of the MNA
        matrix and contains the information regarding which nodes fail to
        converge in the circuit. Ie. ``if convergence_by_node[j] == False``,
        node ``j`` has a convergence problem. This may significantly help
        debugging non-convergent circuits.

    """
    # OLD COMMENT: FIXME REWRITE: solve through newton
    # problem is F(x)= mna*x +H(x) = 0
    # H(x) = N + T(x)
    # lets say: J = dF/dx = mna + dT(x)/dx
    # J*dx = -1*(mna*x+N+T(x))
    # dT/dx e' lo jacobiano -> g_eq (o gm)
    # print_steps = False
    # locked_nodes = get_locked_nodes(element_list)
    mna_size = mna.shape[0]
    nonlinear_circuit = circ.is_nonlinear()
    tick = ticker.ticker(increments_for_step=1)
    tick.display(print_steps)
    if x is None:
        # if no guess was specified, its all zeros
        x = np.zeros((mna_size, 1))
    else:
        if x.shape[0] != mna_size:
            raise ValueError("x0s size is different from expected: got "
                             "%d-elements x0 with an MNA of size %d" %
                             (x.shape[0], mna_size))
    if T is None:
        printing.print_warning(
            "dc_analysis.mdn_solver called with T==None, setting T=0. BUG or no sources in circuit?")
        T = np.zeros((mna_size, 1))

    sparse = mna_size > options.dense_matrix_limit
    # We allocate the matrices once and then reuse them
    if sparse:
        mna = scipy.sparse.coo_matrix(mna)
        J = scipy.sparse.lil_matrix((mna_size, mna_size))
    else:
        J = np.zeros((mna_size, mna_size))
    Tx = np.zeros((mna_size, 1))
    converged = False
    iteration = 0
    while iteration < MAXIT:  # newton iteration counter
        iteration += 1
        tick.step()
        if nonlinear_circuit:
            # build dT(x)/dx (stored in J) and Tx(x)
            J[:, :] = 0.0
            Tx[:, 0] = 0.0
            for elem in circ:
                if elem.is_nonlinear:
                    _update_J_and_Tx(J, Tx, x, elem, time)
        residuo = mna.dot(x) + T + nonlinear_circuit*Tx

        if sparse:
            lu = scipy.sparse.linalg.splu(scipy.sparse.csc_matrix(mna + nonlinear_circuit*J))
            dx = lu.solve(-residuo)
        else:
            dx = np.linalg.solve(mna + nonlinear_circuit*J, -residuo)
        x = x + get_td(dx, locked_nodes, n=iteration) * dx
        if not nonlinear_circuit:
            converged = True
            break
        elif convergence_check(x, dx, residuo, nv - 1)[0]:
            converged = True
            break
        # if vector_norm(dx) == np.nan: #Overflow
        #   raise OverflowError
    tick.hide(print_steps)
    if debug and not converged:
        # re-run the convergence check, only this time get the results
        # by node, so we can show to the users which nodes are misbehaving.
        converged, convergence_by_node = convergence_check(
            x, dx, residuo, nv - 1, debug=True)
    else:
        convergence_by_node = []
    return (x, residuo, converged, iteration, convergence_by_node)


def _update_J_and_Tx(J, Tx, x, elem, time):
    out_ports = elem.get_output_ports()
    for index in range(len(out_ports)):
        n1, n2 = out_ports[index]
        n1m1, n2m1 = n1 - 1, n2 - 1
        dports = elem.get_drive_ports(index)
        v_dports = []
        for port in dports:
            v = 0.  # build v: remember we removed the 0 row and 0 col of mna -> -1
            if port[0]:
                v = v + x[port[0] - 1, 0]
            if port[1]:
                v = v - x[port[1] - 1, 0]
            v_dports.append(v)
        if hasattr(elem, 'gstamp') and hasattr(elem, 'istamp'):
            iis, gs = elem.gstamp(v_dports, time)
            J[iis] += gs.reshape(-1)
            iis, i = elem.istamp(v_dports, time)
            Tx[iis] += i.reshape(-1)
            continue
        if n1 or n2:
            iel = elem.i(index, v_dports, time)
        if n1:
            Tx[n1m1, 0] = Tx[n1m1, 0] + iel
        if n2:
            Tx[n2m1, 0] = Tx[n2m1, 0] - iel
        for iindex in range(len(dports)):
            if n1 or n2:
                g = elem.g(index, v_dports, iindex, time)
            if n1:
                if dports[iindex][0]:
                    J[n1m1, dports[iindex][0] - 1] += g
                if dports[iindex][1]:
                    J[n1m1, dports[iindex][1] - 1] -= g
            if n2:
                if dports[iindex][0]:
                    J[n2m1, dports[iindex][0] - 1] -= g
                if dports[iindex][1]:
                    J[n2m1, dports[iindex][1] - 1] += g


def get_td(dx, locked_nodes, n=-1):
    """Calculates the damping coefficient for the Newthon method.

    The damping coefficient is choosen as the lowest between:

    - the damping required for the first NR iterations, a parameter which is set
      through the integer ``options.nr_damp_first_iters``.
    - If ``options.nl_voltages_lock`` evaluates to ``True``, the biggest damping
      factor that keeps the change in voltage across the locked nodes pairs less
      than the maximum variation allowed, set by:
      ``(options.nl_voltages_lock_factor * Vth)``
    - Unity.

    **Parameters:**

    dx : ndarray
        The undamped increment returned by the NR solver.
    locked_nodes : list
        A vector of tuples of (internal) nodes that are a port of a non-linear
        component.
    n : int, optional
        The NR iteration counter

    .. note::

        If ``n`` is set to ``-1`` (or any negative value), ``td`` is independent
        from the iteration number and ``options.nr_damp_first_iters`` is ignored.

    **Returns:**

    td : float
        The damping coefficient.

    """

    if not options.nr_damp_first_iters or n < 0:
        td = 1
    else:
        if n < 10:
            td = 1e-2
        elif n < 20:
            td = 0.1
        else:
            td = 1
    td_new = 1
    if options.nl_voltages_lock:
        for (n1, n2) in locked_nodes:
            if n1 != 0:
                if n2 != 0:
                    if abs(dx[n1 - 1, 0] - dx[n2 - 1, 0]) > options.nl_voltages_lock_factor * constants.Vth():
                        td_new = (options.nl_voltages_lock_factor * constants.Vth()) / abs(
                            dx[n1 - 1, 0] - dx[n2 - 1, 0])
                else:
                    if abs(dx[n1 - 1, 0]) > options.nl_voltages_lock_factor * constants.Vth():
                        td_new = (options.nl_voltages_lock_factor * constants.Vth()) / abs(
                            dx[n1 - 1, 0])
            else:
                if abs(dx[n2 - 1, 0]) > options.nl_voltages_lock_factor * constants.Vth():
                    td_new = (options.nl_voltages_lock_factor * constants.Vth()) / abs(
                        dx[n2 - 1, 0])
            if td_new < td:
                td = td_new
    return td


def generate_mna_and_N(circ, verbose=3):
    """Generate the full *unreduced* MNA and N matrices required for an MNA analysis

    We wish to solve the linear stationary MNA problem:

    .. math::

        MNA \\cdot x + N = 0

    If ``nv`` is the number of nodes in the circuit, ``MNA`` is a square matrix
    composed by:

    * ``MNA[:nv, :]``, KCLs ordered by node, from node 0 up to node nv.

    In the above submatrix, we have a voltage part: ``MNA[:nv, :nv]``, where
    each term ``MNA[i, j]`` is due to the (trans-)conductances in between the
    nodes and a current part, ``MNA[:nv, nv:]``, where each term is due to a
    current variable introduced by elements whose current flow is not univocally
    defined by the voltage applied to their port(s).

    * ``MNA[nv:, :]`` are the KVL equations introduced by the above terms.

    ``N`` is similarly partitioned, but it is a vector of size ``(nv,)``.

    **Parameters:**

    circ : circuit instance
        The circuit for which the matrices are to be computed.
    verbose : int, optional
        The verbosity, from 0 (silent) to 6 (debug).

    **Returns:**

    MNA, N : ndarrays
        The MNA matrix and constant term vector computed as per above.

    """
    n_of_nodes = circ.get_nodes_number()
    mna = np.zeros((n_of_nodes, n_of_nodes))
    N = np.zeros((n_of_nodes, 1))
    for elem in circ:
        if elem.is_nonlinear:
            continue
        elif isinstance(elem, devices.Resistor):
            mna[elem.n1, elem.n1] = mna[elem.n1, elem.n1] + elem.g
            mna[elem.n1, elem.n2] = mna[elem.n1, elem.n2] - elem.g
            mna[elem.n2, elem.n1] = mna[elem.n2, elem.n1] - elem.g
            mna[elem.n2, elem.n2] = mna[elem.n2, elem.n2] + elem.g
        elif isinstance(elem, devices.Capacitor):
            pass  # In a capacitor I(V) = 0
        elif isinstance(elem, devices.GISource):
            mna[elem.n1, elem.sn1] = mna[elem.n1, elem.sn1] + elem.alpha
            mna[elem.n1, elem.sn2] = mna[elem.n1, elem.sn2] - elem.alpha
            mna[elem.n2, elem.sn1] = mna[elem.n2, elem.sn1] - elem.alpha
            mna[elem.n2, elem.sn2] = mna[elem.n2, elem.sn2] + elem.alpha
        elif isinstance(elem, devices.ISource):
            if not elem.is_timedependent:  # convenzione normale!
                N[elem.n1, 0] = N[elem.n1, 0] + elem.I()
                N[elem.n2, 0] = N[elem.n2, 0] - elem.I()
            else:
                pass  # vengono aggiunti volta per volta
        elif isinstance(elem, devices.InductorCoupling):
            pass
            # this is taken care of within the inductors
        elif circuit.is_elem_voltage_defined(elem):
            pass
            # we'll add its lines afterwards
        elif isinstance(elem, devices.FISource):
            # we add these last, they depend on voltage sources
            # to sense the current
            pass
        else:
            print("dc_analysis.py: BUG - Unknown linear element. Ref. #28934")
    # process vsources
    # i generatori di tensione non sono pilotabili in tensione: g e' infinita
    # for each vsource, introduce a new variable: the current flowing through it.
    # then we introduce a KVL equation to be able to solve the circuit
    for elem in circ:
        if circuit.is_elem_voltage_defined(elem):
            index = mna.shape[0]  # get_matrix_size(mna)[0]
            mna = utilities.expand_matrix(mna, add_a_row=True, add_a_col=True)
            N = utilities.expand_matrix(N, add_a_row=True, add_a_col=False)
            # KCL
            mna[elem.n1, index] = 1.0
            mna[elem.n2, index] = -1.0
            # KVL
            mna[index, elem.n1] = +1.0
            mna[index, elem.n2] = -1.0
            if isinstance(elem, devices.VSource) and not elem.is_timedependent:
                # corretto, se e' def una parte tempo-variabile ci pensa
                # mdn_solver a scegliere quella giusta da usare.
                N[index, 0] = -1.0 * elem.V()
            elif isinstance(elem, devices.VSource) and elem.is_timedependent:
                pass  # taken care step by step
            elif isinstance(elem, devices.EVSource):
                mna[index, elem.sn1] = -1.0 * elem.alpha
                mna[index, elem.sn2] = +1.0 * elem.alpha
            elif isinstance(elem, devices.Inductor):
                # N[index,0] = 0 pass, it's already zero
                pass
            elif isinstance(elem, devices.HVSource):
                index_source = circ.find_vde_index(elem.source_id)
                mna[index, n_of_nodes+index_source] = 1.0 * elem.alpha
            else:
                print("dc_analysis.py: BUG - found an unknown voltage_def elem.")
                print(elem)
                sys.exit(33)

    # iterate again for devices that depend on voltage-defined ones.
    for elem in circ:
        if isinstance(elem, devices.FISource):
            local_i_index = circ.find_vde_index(elem.source_id, verbose=0)
            mna[elem.n1, n_of_nodes + local_i_index] = mna[elem.n1, n_of_nodes + local_i_index] + elem.alpha
            mna[elem.n2, n_of_nodes + local_i_index] = mna[elem.n2, n_of_nodes + local_i_index] - elem.alpha

    # Seems a good place to run some sanity check
    # for the time being we do not halt the execution
    utilities.check_ground_paths(mna, circ, reduced_mna=False, verbose=verbose)

    # all done
    return (mna, N)


def build_x0_from_user_supplied_ic(circ, icdict):
    """Builds a vector of appropriate (reduced!) size from the values supplied
    in ``icdict``.

    Supplying a custom x0 can be useful:
    - To aid convergence in tough circuits,
    - To start a transient simulation from a particular x0.

    **Parameters:**

    circ: circuit instance
        The circuit the :math:`x_0` is being assembled for
    icdict: dict
        ``icdict`` is a a dictionary assembled as follows:
         - to specify a nodal voltage: ``{'V(node)':<voltage value>}``
           Eg. ``{'V(n1)':2.3, 'V(n2)':0.45, ...}``.
           All unspecified voltages default to 0.
         - to specify a branch current: ``'I(<element>)':<current value>}``
           ie. the elements names are sorrounded by ``I(...)``.
           Eg. ``{'I(L1)':1.03e-3, I(V4):2.3e-6, ...}``
           All unspecified currents default to 0.

    Notes: this simulator uses the standard convention.

    **Returns:**

    x0 : ndarray
        The x0 matrix assembled according to ``icdict``.

    :raises ValueError: whenever a malformed ``icdict`` is supplied.
    """
    Vregex = re.compile("V\s*\(\s*([a-z0-9]+)\s*\)", re.IGNORECASE | re.DOTALL)
    Iregex = re.compile("I\s*\(\s*([a-z0-9]+)\s*\)", re.IGNORECASE | re.DOTALL)
    nv = circ.get_nodes_number()  # number of voltage variables
    voltage_defined_elem_names = \
        [elem.part_id.lower() for elem in circ
         if circuit.is_elem_voltage_defined(elem)]
    ni = len(voltage_defined_elem_names)  # number of current variables
    x0 = np.zeros((nv + ni, 1))
    for label in icdict.keys():
        value = icdict[label]
        if Vregex.search(label):
            ext_node = Vregex.findall(label)[0]
            int_node = circ.ext_node_to_int(ext_node)
            x0[int_node, 0] = value
        elif Iregex.search(label):
            element_name = Iregex.findall(label)[0]
            index = voltage_defined_elem_names.index(element_name.lower())
            x0[nv + index, 0] = value
        else:
            raise ValueError("Unrecognized label " + label)
    return x0[1:, :]


def modify_x0_for_ic(circ, x0):
    """Modifies a supplied x0.

    Several circut elements allow the user to set their own Initial
    Conditions (IC) for either voltage or current, depending on what
    is appropriate for the element considered.

    This method, receives a preliminary ``x0`` value, typically computed
    by an OP analysis and goes through the circuit, looking for ICs and
    setting them in ``x0``.

    Notice it is possible to require ICs that are incompatible with each
    other -- for example supplying different ICs to two parallel caps.
    In that case we try to accommodate the user's requirements in a
    non-strict best-effort kind of way: for this reason, whenever
    multiple ICs are specified, it is best to visually inspect ``x0``
    to check that what you would have expected is indeed what you got.

    **Parameters**

    circ : circuit instance
        The circuit in which the ICs are specified.
    x0 : ndarray or results.op_solution
        The initial value to be modified

    **Returns:**

    x0p : ndarray or results.op_solution
        The modified ``x0``. Notice that we return the same
        kind of object as it was supplied. Additionally,
        the ``results.op_solution`` is a *new* *instance*,
        while the ``ndarray`` is simply the original array
        modified.
    """

    if isinstance(x0, results.op_solution):
        x0 = copy.copy(x0.asarray())
        return_obj = True
    else:
        return_obj = False

    nv = circ.get_nodes_number()  # number of voltage variables
    voltage_defined_elements = [
        x for x in circ if circuit.is_elem_voltage_defined(x)]

    # setup voltages this may _not_ work properly
    for elem in circ:
        if isinstance(elem, devices.Capacitor) and elem.ic or \
                isinstance(elem, diode.diode) and elem.ic:
            x0[elem.n1 - 1, 0] = x0[elem.n2 - 1, 0] + elem.ic

    # setup the currents
    for elem in voltage_defined_elements:
        if isinstance(elem, devices.Inductor) and elem.ic:
            x0[nv - 1 + voltage_defined_elements.index(elem), 0] = elem.ic

    if return_obj:
        xnew = results.op_solution(x=x0, \
            error=np.zeros(x0.shape), circ=circ, outfile=None)
        xnew.netlist_file = None
        xnew.netlist_title = "Self-generated OP to be used as tran IC"
    else:
        xnew = x0

    return xnew


