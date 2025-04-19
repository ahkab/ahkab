# -*- coding: iso-8859-1 -*-
# transient.py
# Transient analysis
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

""" This module provides the methods required to perform a transient analysis.

Our problem can be written as:

.. math::

    D \\cdot dx/dt + MNA \\cdot x + T_v(x) + T_t(t) + N = 0

We need:

    1. :math:`MNA`, the static Modified Nodal Analysis matrix,
    2. :math:`N`, constant DC term,
    3. :math:`T_v(x)`, the non-linear DC term
    4. :math:`T_t(t)`, the time variant term, time dependent-sources,
       to be evaluated at each time step,
    5. The dynamic :math:`D` matrix,
    6. a differentiation method to approximate :math:`dx/dt`.

"""

from __future__ import (unicode_literals, absolute_import,
                        division, print_function)
import sys

import numpy as np
import importlib
import importlib.util

from . import dc_analysis
from . import implicit_euler
from . import ticker
from . import options
from . import circuit
from . import printing
from . import utilities
from . import components
from . import results

# differentiation methods, add them here
IMPLICIT_EULER = "IMPLICIT_EULER"
TRAP = "TRAP"
GEAR1 = "GEAR1"
GEAR2 = "GEAR2"
GEAR3 = "GEAR3"
GEAR4 = "GEAR4"
GEAR5 = "GEAR5"
GEAR6 = "GEAR6"

specs = {'tran':{'tokens':({
                          'label':'tstep',
                          'pos':0,
                          'type':float,
                          'needed':True,
                          'dest':'tstep',
                          'default':None
                         },
                         {
                          'label':'tstop',
                          'pos':1,
                          'type':float,
                          'needed':True,
                          'dest':'tstop',
                          'default':None
                         },
                         {
                          'label':'tstart',
                          'pos':None,
                          'type':float,
                          'needed':False,
                          'dest':'tstart',
                          'default':0
                         },
                         {
                          'label':'uic',
                          'pos':2,
                          'type':float,
                          'needed':False,
                          'dest':'uic',
                          'default':0
                         },
                         {
                          'label':'ic_label',
                          'pos':None,
                          'type':str,
                          'needed':False,
                          'dest':'x0',
                          'default':0
                         },
                         {
                          'label':'method',
                          'pos':None,
                          'type':str,
                          'needed':False,
                          'dest':'method',
                          'default':None
                         }
                        )
               }
           }


def transient_analysis(circ, tstart, tstep, tstop, method=options.default_tran_method, use_step_control=True, x0=None,
                       mna=None, N=None, D=None, outfile="stdout", return_req_dict=None, verbose=3):
    """Performs a transient analysis of the circuit described by circ.

    Parameters:
    circ: circuit instance to be simulated.
    tstart: start value. Better leave this to zero.
    tstep: the maximum step to be allowed during simulation or
    tstop: stop value for simulation
    method: differentiation method: 'TRAP' (default) or 'IMPLICIT_EULER' or 'GEARx' with x=1..6
    use_step_control: the LTE will be calculated and the step adjusted. default: True
    x0: the starting point, the solution at t=tstart (defaults to None, will be set to the OP)
    mna, N, D: MNA matrices, defaulting to None, for big circuits, reusing matrices saves time
    outfile: filename, the results will be written to this file. "stdout" means print out.
    return_req_dict:  to be documented
    verbose: verbosity level from 0 (silent) to 6 (very verbose).

    """
    if outfile == "stdout":
        verbose = 0
    _debug = False
    if options.transient_no_step_control:
        use_step_control = False
    if _debug:
        print_step_and_lte = True
    else:
        print_step_and_lte = False

    method = method.upper() if method is not None else options.default_tran_method
    HMAX = tstep

    #check parameters
    if tstart > tstop:
        printing.print_general_error("tstart > tstop")
        sys.exit(1)
    if tstep < 0:
        printing.print_general_error("tstep < 0")
        sys.exit(1)

    if verbose > 4:
        tmpstr = "Vea = %g Ver = %g Iea = %g Ier = %g max_time_iter = %g HMIN = %g" % \
        (options.vea, options.ver, options.iea, options.ier, options.transient_max_time_iter, options.hmin)
        printing.print_info_line((tmpstr, 5), verbose)

    locked_nodes = circ.get_locked_nodes()

    if print_step_and_lte:
        flte = open("step_and_lte.graph", "w")
        flte.write("#T\tStep\tLTE\n")

    printing.print_info_line(("Starting transient analysis: ", 3), verbose)
    printing.print_info_line(("Selected method: %s" % (method,), 3), verbose)
    #It's a good idea to call transient with prebuilt MNA and N matrix
    #the analysis will be slightly faster (long netlists).
    if mna is None or N is None:
        (mna, N) = dc_analysis.generate_mna_and_N(circ, verbose=verbose)
        mna = utilities.remove_row_and_col(mna)
        N = utilities.remove_row(N, rrow=0)
    elif not mna.shape[0] == N.shape[0]:
        printing.print_general_error("mna matrix and N vector have different number of columns.")
        sys.exit(0)
    if D is None:
        # if you do more than one tran analysis, output streams should be changed...
        # this needs to be fixed
        D = generate_D(circ, (mna.shape[0], mna.shape[0]))
        D = utilities.remove_row_and_col(D)

    # setup x0
    if x0 is None:
        printing.print_info_line(("Generating x(t=%g) = 0" % (tstart,), 5), verbose)
        x0 = np.zeros((mna.shape[0], 1))
        opsol =  results.op_solution(x=x0, error=x0, circ=circ, outfile=None)
    else:
        if isinstance(x0, results.op_solution):
            opsol = x0
            x0 = x0.asarray()
        else:
            opsol =  results.op_solution(x=x0, error=np.zeros((mna.shape[0], 1)), circ=circ, outfile=None)
        printing.print_info_line(("Using the supplied op as x(t=%g)." % (tstart,), 5), verbose)

    if verbose > 4:
        print("x0:")
        opsol.print_short()

    # setup the df method
    printing.print_info_line(("Selecting the appropriate DF ("+method+")... ", 5), verbose, print_nl=False)
    if method == IMPLICIT_EULER:
        from . import implicit_euler as df
    elif method == TRAP:
        from . import trap as df
    elif method == GEAR1:
        from . import gear as df
        df.order = 1
    elif method == GEAR2:
        from . import gear as df
        df.order = 2
    elif method == GEAR3:
        from . import gear as df
        df.order = 3
    elif method == GEAR4:
        from . import gear as df
        df.order = 4
    elif method == GEAR5:
        from . import gear as df
        df.order = 5
    elif method == GEAR6:
        from . import gear as df
        df.order = 6
    else:
        df = import_custom_df_module(method, print_out=(outfile != "stdout"))
        # df is none if module is not found

    if df is None:
        sys.exit(23)

    if not df.has_ff() and use_step_control:
        printing.print_warning("The chosen DF does not support step control. Turning off the feature.")
        use_step_control = False
        #use_aposteriori_step_control = False

    printing.print_info_line(("done.", 5), verbose)

    # setup the data buffer
    # if you use the step control, the buffer has to be one point longer.
    # That's because the excess point is used by a FF in the df module to predict the next value.
    printing.print_info_line(("Setting up the buffer... ", 5), verbose, print_nl=False)
    ((max_x, max_dx), (pmax_x, pmax_dx)) = df.get_required_values()
    if max_x is None and max_dx is None:
        printing.print_general_error("df doesn't need any value?")
        sys.exit(1)
    if use_step_control:
        buffer_len = 0
        for mx in (max_x, max_dx, pmax_x, pmax_dx):
            if mx is not None:
                buffer_len = max(buffer_len, mx)
        buffer_len += 1
        thebuffer = dfbuffer(length=buffer_len, width=3)
    else:
        thebuffer = dfbuffer(length=max(max_x, max_dx) + 1, width=3)
    thebuffer.add((tstart, x0, None)) #setup the first values
    printing.print_info_line(("done.", 5), verbose) #FIXME

    #setup the output buffer
    if return_req_dict:
        output_buffer = dfbuffer(length=return_req_dict["points"], width=1)
        output_buffer.add((x0,))
    else:
        output_buffer = None

    # import implicit_euler to be used in the first iterations
    # this is because we don't have any dx when we start, nor any past point value
    if (max_x is not None and max_x > 0) or max_dx is not None:
        from . import implicit_euler
        first_iterations_number = max_x if max_x is not None else 1
        first_iterations_number = max( first_iterations_number, max_dx+1) \
                                  if max_dx is not None else first_iterations_number
    else:
        first_iterations_number = 0

    printing.print_info_line(("MNA (reduced):", 5), verbose)
    printing.print_info_line((mna, 5), verbose)
    printing.print_info_line(("D (reduced):", 5), verbose)
    printing.print_info_line((D, 5), verbose)

    # setup the initial values to start the iteration:
    x = None
    time = tstart
    nv = circ.get_nodes_number()

    Gmin_matrix = dc_analysis.build_gmin_matrix(circ, options.gmin, mna.shape[0], verbose)

    # lo step viene generato automaticamente, ma non superare mai quello fornito.
    if use_step_control:
        #tstep = min((tstop-tstart)/9999.0, HMAX, 100.0 * options.hmin)
        tstep = min((tstop-tstart)/9999.0, HMAX)
    printing.print_info_line(("Initial step: %g"% (tstep,), 5), verbose)

    if max_dx is None:
        max_dx_plus_1 = None
    else:
        max_dx_plus_1 = max_dx +1
    if pmax_dx is None:
        pmax_dx_plus_1 = None
    else:
        pmax_dx_plus_1 = pmax_dx +1

    # setup error vectors
    aerror = np.zeros((x0.shape[0], 1))
    aerror[:nv-1, 0] = options.vea
    aerror[nv-1:, 0] = options.vea
    rerror = np.zeros((x0.shape[0], 1))
    rerror[:nv-1, 0] = options.ver
    rerror[nv-1:, 0] = options.ier

    iter_n = 0  # contatore d'iterazione
    # when to start predicting the next point
    start_pred_iter = max(*[i for i in (0, pmax_x, pmax_dx_plus_1) if i is not None])
    lte = None
    sol = results.tran_solution(circ, tstart, tstop, op=x0, method=method, outfile=outfile)
    printing.print_info_line(("Solving... ", 3), verbose, print_nl=False)
    tick = ticker.ticker(increments_for_step=1)
    tick.display(verbose > 1)
    while time < tstop:
        if iter_n < first_iterations_number:
            x_coeff, const, x_lte_coeff, prediction, pred_lte_coeff = \
            implicit_euler.get_df((thebuffer.get_df_vector()[0],), tstep, \
            predict=(use_step_control and iter_n >= start_pred_iter))
        else:
            x_coeff, const, x_lte_coeff, prediction, pred_lte_coeff = \
                df.get_df(thebuffer.get_df_vector(), tstep,
                          predict=(use_step_control and
                                   iter_n >= start_pred_iter)
                         )

        if options.transient_prediction_as_x0 and use_step_control and prediction is not None:
            x0 = prediction
        elif x is not None:
            x0 = x

        x1, error, solved, n_iter = dc_analysis.dc_solve(
                                                     mna=(mna + np.multiply(x_coeff, D)),
                                                     Ndc=N,  Ntran=np.dot(D, const), circ=circ,
                                                     Gmin=Gmin_matrix, x0=x0,
                                                     time=(time + tstep),
                                                     locked_nodes=locked_nodes,
                                                     MAXIT=options.transient_max_nr_iter,
                                                     verbose=0
                                                     )

        if solved:
            old_step = tstep #we will modify it, if we're using step control otherwise it's the same
            # step control (yeah)
            if use_step_control:
                if x_lte_coeff is not None and pred_lte_coeff is not None and prediction is not None:
                    # this is the Local Truncation Error :)
                    lte = abs((x_lte_coeff / (pred_lte_coeff - x_lte_coeff)) * (prediction - x1))
                    # it should NEVER happen that new_step > 2*tstep, for stability
                    new_step_coeff = 2
                    for index in range(x.shape[0]):
                        if lte[index, 0] != 0:
                            new_value = ((aerror[index, 0] + rerror[index, 0]*abs(x[index, 0])) / lte[index, 0]) \
                            ** (1.0 / (df.order+1))
                            if new_value < new_step_coeff:
                                new_step_coeff = new_value
                            #print new_value
                    new_step = tstep * new_step_coeff
                    if (options.transient_use_aposteriori_step_control and
                        new_step_coeff <
                        options.transient_aposteriori_step_threshold):
                        #don't recalculate a x for a small change
                        tstep = check_step(new_step, time, tstop, HMAX)
                        #print "Apost. (reducing) step = "+str(tstep)
                        continue
                    tstep = check_step(new_step, time, tstop, HMAX) # used in the next iteration
                    #print "Apriori tstep = "+str(tstep)
                else:
                    #print "LTE not calculated."
                    lte = None
            if print_step_and_lte and lte is not None:
                #if you wish to look at the step. We print just a lte
                flte.write(str(time)+"\t"+str(old_step)+"\t"+str(lte.max())+"\n")
            # if we get here, either aposteriori_step_control is
            # disabled, or it's enabled and the error is small
            # enough. Anyway, the result is GOOD, STORE IT.
            time = time + old_step
            x = x1
            iter_n = iter_n + 1
            sol.add_line(time, x)

            dxdt = np.multiply(x_coeff, x) + const
            thebuffer.add((time, x, dxdt))
            if output_buffer is not None:
                output_buffer.add((x, ))
            tick.step()
        else:
            # If we get here, Newton failed to converge. We need to reduce the step...
            if use_step_control:
                tstep = tstep/5.0
                tstep = check_step(tstep, time, tstop, HMAX)
                printing.print_info_line(("At %g s reducing step: %g s (convergence failed)" % (time, tstep), 5), verbose)
            else: #we can't reduce the step
                printing.print_general_error("Can't converge with step "+str(tstep)+".")
                printing.print_general_error("Try setting --t-max-nr to a higher value or set step to a lower one.")
                solved = False
                break
        if options.transient_max_time_iter and iter_n == options.transient_max_time_iter:
            printing.print_general_error("MAX_TIME_ITER exceeded ("+str(options.transient_max_time_iter)+"), iteration halted.")
            solved = False
            break

    if print_step_and_lte:
        flte.close()

    tick.hide(verbose > 1)

    if solved:
        printing.print_info_line(("done.", 3), verbose)
        printing.print_info_line(("Average time step: %g" % ((tstop - tstart)/iter_n,), 3), verbose)

        if output_buffer:
            ret_value = output_buffer.get_as_matrix()
        else:
            ret_value = sol
    else:
        print("failed.")
        ret_value =  None

    return ret_value

def check_step(tstep, time, tstop, HMAX):
    """Checks the step for several common issues and corrects them.

    The following problems are checked:

    - the step must be shorter than ``HMAX``. In the context of a transient
      analysis, that usually is the time step provided by the user,
    - the step must be equal or shorter than the simulation time left (ie stop
      time minus current time),
    - the step must be longer than ``options.hmin``, the minimum allowable time
      step. If the step goes below this value, convergence problems due to
      machine precision will occur. Typically when this happens, we halt the
      simulation.

    **Parameters:**

    tstep : float
        The time step, in second, that needs to be checked.
    time : float
        The current simulation time.
    tstop : float
        The time at which the simulation ends.
    HMAX : float
        The maximum allowable time step.

    **Returns:**

    tstep : float
        The step provided if it passes the tests, a *shortened* step otherwise.

    :raises ValueError: When the step is shorter than ``option.hmin``.

    """
    if tstep > HMAX:
        tstep = HMAX
    if tstop - time < tstep:
        tstep = tstop - time
    elif tstep < options.hmin:
        printing.print_general_error("Step size too small: "+str(tstep))
        raise ValueError("Step size too small")
    return tstep

def generate_D(circ, shape):
    """Generates the D matrix

    For every time t, the D matrix is used (elsewhere) to solve the following system:

    .. math::

        D dx/dt + MNA x + N + T(x) = 0

    It's easy to set up the KCL law for the voltage unknowns, capacitors
    introduce stamps just like resistors do in the MNA and we know that row 1
    refers to node 1, row 2 refers to node 2, and so on

    Inductors generate, together with voltage sources, ccvs, vcvs, a additional
    line in the MNA matrix, and hence in D too.

    The current flowing through the device gets added to the x vector.

    In the case of an inductors, we have:

    .. math::

        V(n_1) - V(n_2) - V_L = 0

    Where:

    .. math::

        V_L = L dI/dt

    That's 0 (zero) in DC analysis, but not in transient analysis, where it
    needs to be differentiated.

    To understand on which line does the inductor's L*dI/dt go, we use the order
    of the elements in `circuit`: first are all voltage lines, then the current
    ones in the same order of the elements that introduce
    them. Therefore, we need to access the circuit (`circ`).

    **Parameters:**

    circ : circuit instance
        The circuit instance for which the :math:`D` matrix is computed.

    shape : tuple of ints
        The shape of the *reduced* :math:`MNA` matrix, D will be of the same
        shape.

    **Returns:**

    D : ndarray
        The *unreduced* D matrix.
    """
    D = np.zeros((shape[0]+1, shape[1]+1))
    nv = circ.get_nodes_number()# - 1
    i_eq = 0 #each time we find a vsource or vcvs or ccvs, we'll add one to this.
    for elem in circ:
        if circuit.is_elem_voltage_defined(elem) and not isinstance(elem, components.Inductor):
            i_eq = i_eq + 1
        elif isinstance(elem, components.Capacitor):
            n1 = elem.n1
            n2 = elem.n2
            D[n1, n1] = D[n1, n1] + elem.value
            D[n1, n2] = D[n1, n2] - elem.value
            D[n2, n2] = D[n2, n2] + elem.value
            D[n2, n1] = D[n2, n1] - elem.value
        elif isinstance(elem, components.Inductor):
            D[ nv + i_eq, nv + i_eq ] = -1 * elem.value
            # Mutual inductors (coupled inductors)
            # need to add a -M dI/dt where I is the current in the OTHER inductor.
            if len(elem.coupling_devices):
                for cd in elem.coupling_devices:
                    # get id+descr of the other inductor (eg. "L32")
                    other_id_wdescr = cd.get_other_inductor(elem.part_id)
                    # find its index to know which column corresponds to its current
                    other_index = circ.find_vde_index(other_id_wdescr, verbose=0)
                    # add the term.
                    D[ nv + i_eq, nv + other_index ] += -1 * cd.M
            # carry on as usual
            i_eq = i_eq + 1

    if options.cmin > 0:
        cmin_mat = np.eye(shape[0]+1-i_eq)
        cmin_mat[0, 1:] = 1
        cmin_mat[1:, 0] = 1
        cmin_mat[0, 0] = cmin_mat.shape[0]-1
        if i_eq:
            D[:-i_eq, :-i_eq] += options.cmin*cmin_mat
        else:
            D += options.cmin*cmin_mat
    return D

class dfbuffer:
    """This is a LIFO buffer with a method to read it all without deleting the elements.

    Newer entries are added on top of the buffer. It checks the size of the
    added elements, to be sure they are of the same size.

    **Parameters:**

    length : int
        The length of the buffer. Samples are added at index ``0``, shifting all
        the previous samples back to higher indices. Samples at an index equal
        to ``length`` (or higher) are discarded without notice.
    width : int
        The width of the buffer, every time :func:`add` is called, it must be to
        add a tuple of the same length as this parameter.
    """
    _the_real_buffer = None
    _length = 0
    _width  = 0

    def __init__(self, length, width):
        self._the_real_buffer = []
        self._length = length
        self._width = width

    def add(self, atuple):
        """Add a new data point to the buffer.

        **Parameters:**

        atuple : tuple of floats
            The data point to be added. Notice that the length of the tuple must
            agree with the width of the buffer.

        :raises ValueError: if the provided tuple and the buffer width do not
        match.

        """
        if not len(atuple) == self._width:
            raise ValueError("Attempted to add a element of wrong size to the" +
                             "LIFO buffer.")
        self._the_real_buffer.insert(0, atuple)
        if len(self._the_real_buffer) > self._length:
            self._the_real_buffer = self._the_real_buffer[:self._length]

    def get_df_vector(self):
        """Read out the contents of the buffer, without any modification

        This method, in the context of a transient analysis, returns a vector
        suitable for a differentiation formula.

        **Returns:**

        vec : list of tuples
            a list of tuples, each tuple being composed of ``width`` floats. In
            the context of a transient analysis, the list (or vector) conforms
            to the specification of the differentiation formulae.
            That is, the simulator stores in the buffer a list similar to::

                [[time(n), x(n), dx(n)], [time(n-1), x(n-1), dx(n-1)], ...]

        """
        return self._the_real_buffer

    def isready(self):
        """This shouldn't be used to determine if the buffer has enough points to
        use the df _if_ you use the step control.
        In that case, it holds even the points required for the FF.
        """
        if len(self._the_real_buffer) == self._length:
            return True
        else:
            return False

    def get_as_matrix(self):
        for vindex in range(self._width):
            for index in range(len(self._the_real_buffer)):
                if index == 0:
                    single_matrix = self._the_real_buffer[index][vindex]
                else:
                    single_matrix = np.concatenate((self._the_real_buffer[index][vindex], single_matrix), axis=0)
            if vindex == 0:
                complete_matrix = single_matrix
            else:
                complete_matrix = np.concatenate((complete_matrix, single_matrix), axis=1)
        return complete_matrix

def import_custom_df_module(method, print_out):
    """Imports a module that implements differentiation formula through imp.load_module
    Parameters:
    method: a string, the name of the df method module
    print_out: print to stdout some verbose messages

    Returns:
    The df module or None if the module is not found.
    """
    try:
        df = importlib.import_module(importlib.util.find_spec(method.lower()))
        if print_out:
            print("Custom df module "+method.lower()+" loaded.")
    except:
        printing.print_general_error("Unrecognized method: "+method.lower()+".")
        df = None

    return df
