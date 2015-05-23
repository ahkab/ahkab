#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-
# ahkab.py
# The frontend of the simulator
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

Introduction
============

This is the core module of the simulator. It provides helper functions to save
you the need to call directly the functions in most submodules.

Do you have a circuit?
======================

To run a simulation, you'll need a circuit first: a circuit can be described
with a simulation deck or with a circuit object.

Define your circuit by means of a Circuit object
------------------------------------------------

In a Python script, describing the circuit through the
:class:`ahkab.circuit.Circuit` interface is a very versatile a choice.

Refer to :class:`ahkab.circuit.Circuit` for a complete description of the
process and the documentation of several helper functions to assist you in this
task.

You may then jump to :ref:`create-simulation-object`.

Define your circuit by means of a netlist file
----------------------------------------------

The circuit description can also be provided as a text file, also known as
netlist deck, for historical reason. This file will also typically include
simulation and post-processing directives, such as plotting.

The netlist should be described according to the rules in
:doc:`help/Netlist-Syntax`.

If you have a netlist (simulation deck) available, you have several possibilities.

The first, assuming your netlist defines some simulation would be to run it:

* you may call ``ahkab`` from the command line. The command line interface is
  described in :doc:`help/Command-Line-Help`.
* you may call :func:`main` directly from Python. Running the simulation through
  :func:`main` function allows to process the result in Python.

Alternatively, you may parse the netlist through
:func:`ahkab.netlist_parser.parse_circuit`, which will return the circuit
instance, all the simulations defined in the deck and all the post-processing
directives as well.

You may now modify the circuit and simulation objects as you please, or create
new ones, as well as run them as described in the :ref:`run-it` section.

.. _create-simulation-object:

How to create a simulation object
=================================

Next, you need to have a simulation object you would like to run.

The following functions are available to quickly create a simulation object:

.. autosummary::
    new_ac
    new_dc
    new_op
    new_pss
    new_pz
    new_symbolic
    new_tran

Click on one of the above hyperlinks to be taken to the corresponding
documentation section.

.. note::

    The functions above allow you to specify an output file. This is due to two
    main reasons:

    * Saving to a file allows you to keep a copy of the simulation results,
      which you can then inspect at a later time.
    * Simulation results may take an uncomfortably large amount of memory. The
      approach we take is that we save everything to file, and only load the
      data to memory when the user actually accesses it.

    In order for the latter to work when no output file is specified, ``ahkab``
    stores the simulation data in a temporary file provided by your OS. When the
    user exits the Python interpreter (or IPython or debugger), the file is
    removed.

.. _run-it:

Run it!
=======

Once you have a circuit and one or more simulations, it's time to run them!

The following methods are available to do so:

.. autosummary::
    run
    queue

The :func:`run` function will return the results in dictionary form.

Extras
======

The core module also contains a few extra methods which were deemed important
enough to be inserted here.

In particular, the :func:`get_op_x0` method allows the user to quickly compute
an operating point to be used to specify the linearization point for a more
complex analysis and :func:`icmodified_x0` allows the user to modify said
operating point to take into account the user-specified initial conditions in
the circuit description.

Lastly, :func:`set_temperature` can be used to quickly set the simulation
temperature.

All methods in alphabetical order
==================================

"""

from __future__ import (unicode_literals, absolute_import,
                        division, print_function)

import atexit
import copy
import os
import sys
import tempfile

from optparse import OptionParser

import numpy as np
import scipy as sp
import sympy
import tabulate
try:
    import matplotlib
    plotting_available = True
except ImportError:
    plotting_available = False

# analyses
from . import dc_analysis
from . import transient
from . import ac
from . import pss
from . import symbolic
from . import pz

# parser
from . import netlist_parser

# misc
from . import options
from . import constants
from . import utilities

# post-processing
from . import fourier

# data display
from . import plotting
from . import printing

#py3 compat
from .py3compat import text_type

from .__version__ import __version__

global _queue, _x0s, _print

_queue = []
_print = False
_x0s = {None: None}

def new_op(guess=None, x0=None, outfile=None, verbose=0):
    """Assembles an OP analysis and returns the analysis object.

    The analysis itself can then be run with: ``ahkab.run(...)``
    or queued with ``ahkab.queue(...)`` and then run subsequently.

    **Parameters:**

    guess : boolean, optional
        if set to True, the analysis will start from an initial guess,
        hopefully speeding up the convergence of stiff circuits.

    x0 : matrix, optional
        In alternative to the ``guess`` option above, one can provide
        an explicit starting point to the OP algorithm, setting x0 to an opportunely sized
        ``numpy`` array. FIXME mention help method here
        If both x0 and guess are set, x0 takes the precedence.

    outfile : string, optional
        the filename of the output file where the results will be written.
        ``.opinfo`` is automatically added at the end to prevent different
        analyses from overwriting each-other's results.
        If unset or set to None, defaults to ``stdout``, if the simulator was
        called from the command line, otherwise, if the simulator is run from an
        interactive session, a temporary file will be used to store the data.

    verbose : int, optional
        the verbosity level, from 0 (silent, default) to 6 (debug).

    **Returns:**

    an : dict
        the analysis description

    .. seealso:: :func:`run`, :func:`queue`
    """
    if guess is None:
        guess = options.dc_use_guess
    if outfile is None or outfile == 'stdout':
        if options.cli:
            outfile = 'stdout'
        else:
            tmpfile = tempfile.NamedTemporaryFile(suffix='.op', delete=False)
            outfile = tmpfile.name
            tmpfile.close()
            atexit.register(os.remove, outfile)
    else:
        outfile += '.op'
    return {'type': 'op', 'guess': guess, 'x0': x0, 'outfile': outfile, 'verbose': verbose}


def new_dc(start, stop, points, source, sweep_type='LINEAR', guess=True, x0=None,
        outfile=None, verbose=0):
    """Assembles a DC sweep analysis and returns the analysis object.

    The analysis itself can be run with: ``ahkab.run(...)``
    or queued and then run subsequently.

    **Parameters:**

    start : float
        the start value for the sweep.

    stop : float
        the stop value for the sweep (included in the sweep points).

    points : int
        the number of sweep points.

    source : string
        the ``part_id`` of the independent current or voltage source to be swept.

    sweep_type : string, optional
        can be set to either ``options.dc_lin_step`` (linear stepping) or
        ``options.dc_log_step`` (log10 stepping). Defaults to linear.

    guess : boolean, optional
        if set to ``True``, the analysis will start from an initial guess,
        hopefully speeding up the convergence of particularly stiff circuits.

    x0 : ``numpy`` array, optional
        if the ``guess`` option above is not used, one can provide
        a starting point directly, setting ``x0`` to an opportunely sized
        ``numpy`` array.
        If both ``x0`` and ``guess`` are set, ``x0`` takes the precedence.

    outfile : string, optional
        the filename of the output file where the results will be written.
        '.dc' is automatically added at the end to prevent different
        analyses from overwriting each-other's results.
        If unset or set to ``None``, defaults to ``stdout``, if the simulator
        was called from the command line, otherwise, if the simulator is run
        from an interactive session, a temporary file will be used to store the
        data.

    verbose : int, optional
        the verbosity level, from 0 (silent, default) to 6 (debug).

    **Returns:**

    an : dict
        the analysis description

    .. seealso:: :func:`run`, :func:`queue`
    """
    if outfile is None or outfile == 'stdout':
        if options.cli:
            outfile = 'stdout'
        else:
            tmpfile = tempfile.NamedTemporaryFile(suffix='.dc', delete=False)
            outfile = tmpfile.name
            tmpfile.close()
            atexit.register(os.remove, outfile)
    else:
        outfile += '.dc'
    return {
        'type': 'dc', 'start': float(start), 'stop': float(stop), 'step': float(stop - start) / float(points - 1),
        'source': source, 'x0': x0, 'outfile': outfile, 'guess': guess, 'sweep_type': sweep_type,
        'verbose': verbose}


def new_tran(tstart, tstop, tstep, x0='op', method=transient.TRAP,
        use_step_control=True, outfile=None, verbose=0):
    """Assembles a TRAN analysis and returns the analysis object.

    The analysis itself can be run with ``ahkab.run(...)``
    or queued with ``ahkab.queue(...)`` and then run subsequently.

    **Parameters:**

    tstart : float
        the start time for the transient analysis.

    tstop : float
        the stop time.

    tstep :float
        the time step. If the step control is active, this is the
        minimum time step value that will be allowed during simulation.

    x0 : ``numpy`` array, optional
        the optional initial conditions point, :math:`x0 = x(t=0)`.

    method : string , optional
        the differentiation method to be used. Can be set to
        'IMPLICIT_EULER', 'TRAP', 'GEAR4', 'GEAR5' or 'GEAR6'.
        It defaults to 'TRAP'.

    use_step_control : boolean, optional
        Whether ste control should be enabled or not. if set to ``False``, the
        differentiation method will use a fixed time step equal to ``tstep``.

    outfile : string, optional
        the filename of the output file where the results will be written.
        '.tran' is automatically added at the end to prevent different
        analyses from overwriting each-other's results.
        If unset or set to ``None``, defaults to ``stdout``, if the simulator
        was called from the command line, otherwise, if the simulator is run
        from an interactive session, a temporary file will be used to store the
        data.

    verbose : int, optional
        the verbosity level, from 0 (silent, default) to 6 (debug).

    **Returns:**

    an : dict
        the analysis description

    .. seealso:: :func:`run`, :func:`queue`
    """
    if outfile is None or outfile == 'stdout':
        if options.cli:
            outfile = 'stdout'
        else:
            tmpfile = tempfile.NamedTemporaryFile(suffix='.tran', delete=False)
            outfile = tmpfile.name
            tmpfile.close()
            atexit.register(os.remove, outfile)
    else:
        outfile += '.tran'
    return {"type": "tran", "tstart": tstart, "tstop": tstop, "tstep": tstep,
            "method": method, "use_step_control": use_step_control, 'x0': x0,
            'outfile': outfile, 'verbose': verbose}


def new_ac(start, stop, points, x0='op', sweep_type='LOG', outfile=None, verbose=0):
    """Assembles an AC analysis and returns the analysis object.

    The analysis itself can be run with ``ahkab.run(...)``
    or queued with ``ahkab.queue(...)`` and then run subsequently.

    **Parameters:**

    start : float
        the start angular frequency, :math:`\omega _{start}`.

    stop : float
        the stop angular frequency, :math:`\omega _{stop}` (included in the
        sweep).

    points : float
        the number of points to be used the discretize the
        `[start, stop]` interval.

    x0 : string or ndarray, optional
        The linearization point for the AC analysis. If set to 'op' (default),
        the latest Operating point analysis will be used. Otherwise, you may
        supply your own linearization point in ndarray format.

    sweep_type : string, optional
        It can be set to either ``options.ac_lin_step`` (linear stepping) or
        ``options.ac_log_step`` (log10 stepping). Defaults to logarithmic
        stepping.

    outfile : string, optional
        the filename of the output file where the results will be written.
        '.ac' is automatically added at the end to prevent different
        analyses from overwriting each-other's results.
        If unset or set to ``None``, defaults to ``stdout``, if the simulator
        was called from the command line, otherwise, if the simulator is run
        from an interactive session, a temporary file will be used to store the
        data.

    verbose : int, optional
        the verbosity level, from 0 (silent, default) to 6 (debug).

    **Returns:**

    an : dict
        the analysis object (a dict)

    .. seealso:: :func:`run`, :func:`queue`
    """
    if outfile is None or outfile == 'stdout':
        if options.cli:
            outfile = 'stdout'
        else:
            tmpfile = tempfile.NamedTemporaryFile(suffix='.ac', delete=False)
            outfile = tmpfile.name
            tmpfile.close()
            atexit.register(os.remove, outfile)
    else:
        outfile += '.ac'
    return {
        'type': 'ac', 'start': start, 'stop': stop, 'points': points,
        'sweep_type': sweep_type, 'x0': x0, 'outfile': outfile,
        'verbose': verbose}


def new_pss(period, x0=None, points=None, method=options.BFPSS, autonomous=False,
            outfile=None, verbose=0):
    """Assembles a Periodic Steady State (PSS) analysis and returns the analysis object.

    The analysis itself can be run with: ``ahkab.run(...)`` or queued with
    ``ahkab.queue(...)`` and then run subsequently.

    **Parameters:**

    period : float
        the time period of the solution, in seconds. This value is required,
        autonomous circuits are currently unsupported.

    x0 : ``numpy`` array, optional
        the starting point solution, used at :math:`t=0`.

    points : int, optional
        the number of points to use to discretize the PSS solution. If not set,
        if method is 'shooting', defaults to ``options.shooting_default_points``

    method : string, optional
        The method to be employed to attempt a PSS solution of the circuit.
        It can be either ``ahkab.BFPSS`` or ``ahkab.SHOOTING``.

    autonomous : bool, optional
        Whether the circuit is autonomous or not.
        Non-autonomous circuits are currently unsupported!

    mna, Tf, D : ``numpy`` arrays, optional
        The matrices to be used to solve the circuit.
        They are optional, if they have already been computed, reusing them saves time.

    outfile : string, optional
        The filename of the output file where the results will be written.
        '.tran' is automatically added at the end to prevent different
        analyses from overwriting each-other's results.
        If unset or set to ``None``, defaults to ``stdout``, if the simulator
        was called from the command line, otherwise, if the simulator is run
        from an interactive session, a temporary file will be used to store the
        data.

    verbose : int, optional
        The verbosity level, from 0 (silent, default) to 6 (debug).

    **Returns:**

    an : dict
        the analysis object (a dict)

    .. seealso:: :func:`run`, :func:`queue`
    """
    if outfile is None or outfile == 'stdout':
        if options.cli:
            outfile = 'stdout'
        else:
            tmpfile = tempfile.NamedTemporaryFile(suffix='.'+method.lower(), delete=False)
            outfile = tmpfile.name
            tmpfile.close()
            atexit.register(os.remove, outfile)
    else:
        outfile += '.' + method.lower()
    return {
        'type': "pss", "method": method, 'period': period, 'points': points,
        'autonomous': autonomous, 'x0': x0, 'outfile': outfile, 'verbose': verbose}


def new_pz(input_source=None, output_port=None, shift=0.0, MNA=None, outfile=None,
           x0='op', verbose=0):
    """Assembles a Pole-Zero analysis and returns the analysis object.

    The analysis itself can be run with: ``ahkab.run(...)`` or queued with
    ``ahkab.queue(...)`` and then run subsequently.

    **Parameters:**

    input_source : str or instance
        the input source for zero calculation

    output_port : tuple or single node
        the output port. If it is composed of only one node, then the
        second node is assumed to be GND.

    shift : float, optional
        Perform the calculation at a shifted freq ``shift``.

    MNA : ndarray, optional
        the numpy matrix to be used to solve the circuit.
        It is optional, but, if it's already been computed,
        reusing it will save time.

    outfile : string, optional
        The filename of the output file where the results will be written.
        '.pz' is automatically added at the end to prevent different
        analyses from overwriting each-other's results.
        If unset or set to ``None``, defaults to ``stdout``, if the simulator
        was called from the command line, otherwise, if the simulator is run
        from an interactive session, a temporary file will be used to store the
        data.

    x0 : ``numpy`` array or str, optional
        the optional linearization point. If set to a string, it must be
        the result of an .OP analysis (use ``'op'``) or an .IC condition
        defined in the netlist. It has no effect on linear circuits.

    verbose : int, optional
        The verbosity level, from 0 (silent, default) to 6 (debug).

    **Returns:**

    an : the analysis description object, a dict instance.
    """
    if outfile is None or outfile == 'stdout':
        if options.cli:
            outfile = 'stdout'
        else:
            tmpfile = tempfile.NamedTemporaryFile(suffix='.pz', delete=False)
            outfile = tmpfile.name
            tmpfile.close()
            atexit.register(os.remove, outfile)
    else:
        outfile += '.pz'
    return {'type': "pz", 'input_source':input_source, 'x0':x0,
            'output_port':output_port, 'shift':shift, 'MNA':MNA,
            'outfile':outfile, 'verbose':verbose}


def new_symbolic(source=None, ac_enable=True, r0s=False, subs=None, outfile=None, verbose=0):
    """Assembles a Symbolic analysis and returns the analysis object.

    The analysis itself can be run with ``ahkab.run(...)``
    or queued with ``ahkab.queue(...)`` and then run subsequently.

    **Parameters:**

    source : string, optional
        if ``source`` is set, the transfer function between the current or
        voltage source ``source`` and each circuit unknown will be
        evaluated, with symbolic evaluation of DC gain, poles and zeros.
        ``source`` is to be set to the ``part_id`` of an independent current
        or voltage source in the circuit, eg. ``'V1'`` or ``'Iin'``. This
        computation should be avoided for large circuit, as indiscriminate
        transfer function, gain and singularities evaluation in large
        circuits can result in very long run times and needs a significant
        amount of RAM, on top of an already resource intensive symbolic
        analysis.  We suggest manually evaluating selected transfer
        functions of interest instead.

    ac_enable : bool, optional
        If set to ``True`` (default), the frequency-dependent elements will
        be considered, otherwise the algorithm will focus on
        low frequency solutions, where all capacitors are replaced with open
        circuits and all inductors are short circuits, usually providing a much
        easier circuit.

    r0s : bool, optional
        If set to ``True``, the finite output conductances of transistors
        ``go`` (where :math:`go = 1/r_0`) will be taken into account,
        otherwise they will be considered infinite (default).
        The finite output conductances generally introduce a significant
        additional complexity in large circuits, sometimes of interest to the
        designer, sometimes simply introducing 2nd and 3rd order effects of
        little-to-no interest, which would produce no significant contribution
        in a numerical analysis, but come at a high computation price in a
        symbolic analysis.
        A possible approach in those cases may be disabling this option and
        explicitly introducing additional conductances where deemed of interest.

    subs : dict, optional
        ``subs`` is a dictionary of substitutions to be performed before
        attempting to solve the circuit. For example, if two
        resistances ``R1`` and ``R2`` are to be equal, set ``subs={'R2':'R1'}``
        and ``R1`` will be replaced by an instance of ``R2``. This may
        simplify the solution (or allow finding one in reasonable
        time for complex circuits).

    outfile : string, optional
        The filename of the output file where the results will be written.
        '.symbolic' is automatically added at the end to prevent different
        analyses from overwriting each-other's results.
        If unset or set to ``None``, defaults to ``stdout``, if the simulator
        was called from the command line, otherwise, if the simulator is run
        from an interactive session, a temporary file will be used to store the
        data.

    verbose : int, optional
        The verbosity level, from 0 (silent, default) to 6 (debug).

    **Returns:**

    an : dict
        the analysis description

    .. seealso:: :func:`run`, :func:`queue`
    """
    if outfile is None or outfile == 'stdout':
        if options.cli:
            outfile = 'stdout'
        else:
            tmpfile = tempfile.NamedTemporaryFile(suffix='.symbolic', delete=False)
            outfile = tmpfile.name
            tmpfile.close()
            atexit.register(os.remove, outfile)
    else:
        outfile += '.symbolic'
    return {
        'type': "symbolic", 'source': source, 'ac_enable': ac_enable, 'r0s': r0s, 'subs': subs,
        'outfile': outfile, 'verbose': verbose}


def queue(*analysis):
    """Queue one or more analyses to execute them subsequently with :func:`run`.

    **Parameters**

    analysis : one or more analysis descriptions.
        The analyses to be queued.

    **Returns:**

    None
    """
    global _queue
    for an in analysis:  # let's hope the user knows what he's doing!
        _queue += [copy.deepcopy(an)]


def run(circ, an_list=None):
    """Run analyses on a circuit.

    **Parameters:**

    circ : circuit instance
        The circuit to be simulated.

    an_queue : list, optional
        the list of analyses to be performed. If unset, it defaults to those
        queued with ``queue``.

    **Returns:**

    results : dict
        the results (in dict form)

    .. seealso:: :func:`queue`
    """
    results = {}

    if not an_list:
        an_list = _queue
    else:
        an_list = copy.deepcopy(an_list)
        if type(an_list) == tuple:
            an_list = list(an_list)
        elif type(an_list) == dict:
            an_list = [an_list] # run(mycircuit, op1)

    while len(an_list):
        an_item = an_list.pop(0)
        an_type = an_item.pop('type')
        if 'x0' in an_item and isinstance(an_item['x0'], text_type):
            printing.print_warning("%s has x0 set to %s, unavailable. Using 'None'." %
                                   (an_type.upper(), an_item['x0']))
            an_item['x0'] = None
        r = analysis[an_type](circ, **an_item)
        results.update({an_type: r})
        if an_type == 'op':
            _x0s.update({'op': r})
            _x0s.update({'op+ic': icmodified_x0(circ, r)})
            _handle_netlist_ics(circ, an_list, ic_list=[])
    return results


def new_x0(circ, icdict):
    """Builds an ``x0`` matrix from user supplied values.

    Supplying a custom x0 can be useful:
    - To aid convergence in tough circuits.
    - To start a transient simulation from a particular x0

    **Parameters:**

    circ : circuit instance
        The circuit
    icdict : dict
        a dictionary specifying the node voltages and branch currents,
        where appropriate, in V and A, respectively, assembled as shown
        in the next section.
        All unspecified node voltages default to ``0`` V and all
        unspecified currents default to 0.

    The user-specified values are to be provided as follows:

    - to specify a nodal voltage: ``{'V(node)':<voltage value>}``
    - to specify a branch current: ``'I(<element>)':<current value>}``


    Examples:

    - ``{'V(n1)':2.3, 'V(n2)':0.45, ...}``
    - ``{'I(L1)':1.03e-3, I(V4):2.3e-6, ...}``

    .. note::
        This simulator uses the normal convention, also known as the
        `Passive sign convention
        <https://en.wikipedia.org/wiki/Passive_sign_convention>`_.

    **Returns:**

    x0 : ``numpy`` array
        The assembled x0.
    """

    return dc_analysis.build_x0_from_user_supplied_ic(circ, icdict)


def icmodified_x0(circ, x0):
    """Modify ``x0`` to take into account the ICs in the circuit.

    **Parameters:**

    circ : circuit instance
        The circuit instance from which the initial conditions are to be
        extracted.

    x0 : ``numpy`` array
        The vector to which the initial conditions are to be applied.
    """
    return dc_analysis.modify_x0_for_ic(circ, x0)


def get_op_x0(circ):
    """Shorthand to specify and run an OP analysis to get a linearization
    point.

    **Parameters:**

    circ : circuit instance
        The circuit instance for which the linearization point is sought.

    **Returns:**

    x0 : an OP solution object
        The linearization point.
    """
    return run(circ, [new_op()])


def set_temperature(T):
    """Set the simulation temperature, in Celsius."""
    T = float(T)
    if T > 300:
        printing.print_warning("The temperature will be set to %f \xB0 C.")
    constants.T = utilities.Celsius2Kelvin(T)


def process_postproc(postproc_list, title, results, outfilename):
    """Runs the post-processing operations, such as plotting.

    Not meant for end users.

    deprecated in 0.10

    **Parameters:**

    postproc_list : list,
        list of post processing operations

    title : string
        the deck title

    results: dict
        the results to be plotted (which may include including ones that are not needed too).

    outfilename: string
        if the plots are saved to disk, this is the filename without extension
    """
    show_plots = False
    skip_plots = False
    for postproc in postproc_list:
        if postproc['type'] == "plot" and outfilename == 'stdout':
            skip_plots = True
            printing.print_warning("Plotting and printing the results to " +
                                   "stdout are incompatible options. " +
                                   "Plotting skipped.")
            break
    for index, postproc in enumerate(postproc_list):
        if postproc['type'] == "plot" and not skip_plots:
            plotting.plot_results(title, postproc["l2l1"],
                                  results[postproc["analysis"]],
                                  "%s-%d.%s" % (outfilename, index,
                                                options.plotting_outtype))
            show_plots = True
        elif postproc['type'] == "four":
            for v in postproc['variables']:
                f, F, THD = fourier.fourier(v, results['tran'],
                                            postproc['fund'])
                printing.print_fourier(v, f, F, THD, outfile='stdout')
        elif postproc['type'] == "fft":
            postproc.pop('type') # get rid of spurious entry
            f, F, THD = fourier.spicefft(tran_results=results['tran'],
                                         **postproc)
            printing.print_spicefft(postproc['label'], f, F, THD,
                                    uformat=postproc['uformat'],
                                    window=postproc['window'],
                                    outfile=outfilename +
                                    '.lis'*(outfilename!='stdout'))
    if show_plots and options.plotting_show_plots:
        plotting.show_plots()

analysis = {'op': dc_analysis.op_analysis, 'dc': dc_analysis.dc_analysis,
            'tran': transient.transient_analysis, 'ac': ac.ac_analysis,
            'pss': pss.pss_analysis, 'symbolic': symbolic.symbolic_analysis,
            'temp': set_temperature, 'pz':pz.calculate_singularities}


def main(filename, outfile="stdout", verbose=3):
    """Method to call ahkab from a Python script with a netlist file.

    **Parameters:**

    filename : string
        The netlist filename.

    outfile : string, optional
        The outfiles base name, the suffixes shown below will be added.
        With the exception of the magic value ``stdout`` which causes
        ahkab to print out instead of to disk.

    verbose : int, optional
        the verbosity level, from 0 (silent) to 6 (debug).
        It defaults to 3, the same as running ahkab through its command
        line interface.

    Filename suffixes, for each analysis:

    - Alternate Current (AC): ``.ac``
    - Direct Current (DC): ``.dc``
    - Operating Point (OP): ``.opinfo``
    - Periodic Steady State (PSS): ``.pss``
    - Pole-zero Analysis (PZ): ``.pz``
    - TRANsient (TRAN): ``.tran``
    - Symbolic: ``.symbolic``

    **Returns:**

    res : dict
        A dictionary containing the computed results.
    """
    printing.print_info_line(
        ("This is ahkab %s running with:" % (__version__), 6), verbose)
    printing.print_info_line(
        ("  Python %s" % (sys.version.split('\n')[0],), 6), verbose)
    printing.print_info_line(("  Numpy %s" % (np.__version__), 6), verbose)
    printing.print_info_line(("  Scipy %s" % (sp.__version__), 6), verbose)
    printing.print_info_line(("  Sympy %s" % (sympy.__version__), 6), verbose)
    printing.print_info_line(("  Tabulate %s" % (tabulate.__version__), 6), verbose)
    if plotting_available:
        printing.print_info_line(("  Matplotlib %s" % (matplotlib.__version__),
                                  6), verbose)
        printing.print_info_line(("  -> backend: %s" %
                                 (matplotlib.get_backend()), 6), verbose)
        printing.print_info_line(("  -> matplotlibrc: %s" %
                                 (matplotlib.matplotlib_fname()), 6), verbose)
    else:
        printing.print_info_line(
            ("  Matplotlib not found.", 6), verbose)


    read_netlist_from_stdin = (filename is None or filename == "-")
    (circ, directives, postproc_direct) = netlist_parser.parse_circuit(
        filename, read_netlist_from_stdin)

    printing.print_info_line(("Checking circuit for common mistakes...", 6),
                             verbose, print_nl=False)
    check, reason = utilities.check_circuit(circ)
    if not check:
        printing.print_general_error(reason)
        sys.exit(3)
    printing.print_info_line(("done.", 6), verbose)

    if verbose > 3 or _print:
        print("Parsed circuit:")
        print(circ)
        print("Models:")
        for m in circ.models:
            circ.models[m].print_model()
        print("")

    ic_list = netlist_parser.parse_ics(directives)
    _handle_netlist_ics(circ, an_list=[], ic_list=ic_list)
    results = {}
    for an in netlist_parser.parse_analysis(circ, directives):
        if 'outfile' not in list(an.keys()) or not an['outfile']:
            an.update(
                {'outfile': outfile + ("." + an['type']) * (outfile != 'stdout')})
        if 'verbose' in list(an.keys()) and (an['verbose'] is None or an['verbose'] < verbose) \
           or not 'verbose' in list(an.keys()):
            an.update({'verbose': verbose})
        _handle_netlist_ics(circ, [an], ic_list=[])
        if verbose >= 4:
            printing.print_info_line(("Requested an.:", 4), verbose)
            printing.print_analysis(an)
        results.update(run(circ, [an]))

    postproc_list = netlist_parser.parse_postproc(circ, postproc_direct)
    if len(postproc_list) > 0 and len(results):
        process_postproc(postproc_list, circ.title, results, outfile)

    return results


def _handle_netlist_ics(circ, an_list, ic_list):
    for ic in ic_list:
        ic_label = list(ic.keys())[0]
        icdict = ic[ic_label]
        _x0s.update({ic_label: new_x0(circ, icdict)})
    for an in an_list:
        if 'x0' in an and isinstance(an['x0'], text_type):
            if an['x0'] in list(_x0s.keys()):
                an['x0'] = _x0s[an['x0']]
            elif an_list.index(an) == 0:
                raise ValueError(("The x0 '%s' is not available." % an["x0"]) +\
                                 (an['x0'] == 'op' or an['x0'] == 'op+ic')*
                                 " Perhaps you forgot to define an .OP?")

