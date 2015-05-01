# -*- coding: iso-8859-1 -*-
# printing.py
# Printing module
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
This is the printing module of the simulator.
Using its functions, the output will be somewhat uniform.

The functions defined in this module can be divided in the following groups:

- :ref:`info-functions`: functions to print information, errors and warnings to
  the user during command-line execution.
- :ref:`netlist_syntax_printing`: functions to print conformingly to the netlist
  syntax, often to show information to the user for debugging purposes.
- :ref:`convenience_functions`: functions to abstract low level issues such as
  Unicode handling of text and number printing formats,
- :ref:`tabular_functions`: functions to format and display data into tables,
  which we provide straight from the ``tabulate`` module.
- :ref:`printing_analysis_results` in a consistent fashion.


.. _info-functions:

Informative functions
=====================

.. autosummary::
    print_general_error
    print_info_line
    print_parse_error
    print_warning

.. _netlist_syntax_printing:

Printing netlist lines
======================

.. autosummary::
    print_analysis

.. _printing_analysis_results:

Printing analysis results
=========================

.. autosummary::
    print_symbolic_equations
    print_symbolic_results
    print_symbolic_transfer_functions

.. _convenience_functions:

Convenience functions
---------------------

.. autosummary::
    open_utf8
    printoptions

.. _tabular_functions:

Tabular formatting of data
--------------------------

.. autosummary::
    table

All functions in alphabetical order
===================================

"""

from __future__ import (unicode_literals, absolute_import,
                        division, print_function)

import contextlib
import sys
import os

import tabulate as _tabulate
import numpy as np

from . import options
from . import py3compat

if py3compat.PY2:
    import codecs
    if not py3compat.IPYTHON:
        UTF8Writer = codecs.getwriter('utf8')
        sys.stdout = UTF8Writer(sys.stdout)

def open_utf8(filename):
    """Get a file handle wrapped in a UTF-8 writer

    The file is opened in ``w`` mode.

    **Parameters:**

    filename : string
        The file name, just like you would pass to Python's built-in ``open()``
        method.

    **Returns:**

    fp : codecs.UTF8Writer object
        The wrapped file pointer.
    """
    fp = open(filename, 'w')
    if py3compat.PY2:
        UTF8Writer = codecs.getwriter('utf8')
        fp = UTF8Writer(fp)
    return fp

def print_analysis(an):
    """Prints an analysis to ``stdout`` in the netlist syntax

    **Parameters:**

    an : dict
        An analysis description in dictionary format.

    """
    if an["type"] == "op":
        print(".op", end="")
        for x in an:
            if x == 'type' or x == 'outfile' or x == 'verbose':
                continue
            print(" %s=%s" % (x, an[x]), end="")
        print("")
    elif an["type"] == "dc":
        print(".dc %(source)s start=%(start)g stop=%(stop)g step=%(step)g type=%(sweep_type)s" % an)
    elif an["type"] == "tran":
        sys.stdout.write(".tran tstep=" + str(an["tstep"]) + " tstop=" + str(
            an["tstop"]) + " tstart=" + str(an["tstart"]))
        if an["method"] is not None:
            print(" method=" + an["method"])
        else:
            print("")
    elif an["type"] == "shooting":
        sys.stdout.write(".shooting period=" + str(
            an["period"]) + " method=" + str(an["method"]))
        if an["points"] is not None:
            sys.stdout.write(" points=" + str(an["points"]))
        if an["step"] is not None:
            sys.stdout.write(" step=" + str(an["step"]))
        print(" autonomous=", an["autonomous"])


def print_general_error(description, print_to_stdout=False):
    """Prints an error message to ``stderr``

    **Parameters:**

    description : str
        The error description.

    print_to_stdout : bool, optional
        When set to ``True``, printing to ``stdout`` instead of ``stderr``.
        Defaults to ``False``.

    """
    the_error_message = "E: " + description
    if print_to_stdout:
        print(the_error_message)
    else:
        sys.stderr.write(the_error_message + "\n")
    return None


def print_warning(description, print_to_stdout=False):
    """Prints a warning message to ``stderr``

    **Parameters:**

    description : str
        The warning message.

    print_to_stdout : bool, optional
        When set to ``True``, printing to ``stdout`` instead of ``stderr``.
        Defaults to ``False``.

    """
    the_warning_message = "W: " + description
    if print_to_stdout:
        print(the_warning_message)
    else:
        sys.stderr.write(the_warning_message + "\n")
    return None


def print_info_line(msg_relevance_tuple, verbose, print_nl=True):
    """Conditionally print out a message

    **Parameters:**

    msg_relevance_tuple : sequence
        A tuple or list made of ``msg`` and ``importance``, where ``msg`` is a
        string, containing the information to be displayed to the user, and
        ``importance``, an integer, is its importance level. Zero corresponds to
        the highest possible importance level, which is always printed out by
        the simple algorithm discussed below.
    verbose : int
        The verbosity level of the program execution. Admissible levels are in
        the 0-6 range.
    print_nl : boolean, optional
        Whether a new line character should be appended or not to the string
        ``msg`` described above, if it's printed out. Defaults to ``True``.

    **Algorithm selecting when to print:**

    The message ``msg`` is printed out if the verbosity level is greater or
    equal than its importance.
    """
    msg, relevance = msg_relevance_tuple
    if verbose >= relevance:
        with printoptions(precision=options.print_precision,
                          suppress=options.print_suppress):
            if print_nl:
                print(msg)
            else:
                print(msg, end=' ')
                sys.stdout.flush()
    # else: suppressed.


def print_parse_error(nline, line, print_to_stdout=False):
    """Prints a parsing error to ``stderr``

    **Parameters:**

    nline : int
        The number of the line on which the error occurred.

    line : str
        The line of the file with the error.

    print_to_stdout : bool, optional
        When set to ``True``, printing to ``stdout`` instead of ``stderr``.
        Defaults to ``False``.

    """
    print_general_error(
        "Parse error on line " + str(nline) + ":", print_to_stdout)
    if print_to_stdout:
        print(line)
    else:
        sys.stderr.write(line + "\n")
    return None


def print_symbolic_results(x):
    """Print out symbolic results

    **Parameters:**

    x : dict
        A dictionary composed of elements like ``{v:expr}``,
        where ``v`` is a circuit variable and ``expr`` is the ``sympy``
        expression corresponding to it, as found by the solver.

    """
    keys = list(x.keys())
    keys.sort(lambda x, y: cmp(str(x), str(y)))
    for key in keys:
        print(str(key) + "\t = " + str(x[key]))


def print_symbolic_transfer_functions(x):
    """Print symbolic transfer functions

    **Parameters:**

    x : dict
        A dictionary of dictionaries. Each top level dictionary is a
        symbol : symbolic transfer function pair, eg. ``{vo/vin:<tf>}``.
        Each transfer function (``<tf>``) is itself a dictionary, having as keys
        the following strings: ``'gain'``, corresponding to the complete
        symbolic TF expression, ``'gain0'``, corresponding to the DC gain and
        ``'poles'`` and ``'zeros'``, corresponding to lists of symbolic
        expressions of the singularities.
    """
    keys = list(x.keys())
    keys.sort(lambda x, y: cmp(str(x), str(y)))
    for key in keys:
        print(str(key) + " = " + str(x[key]['gain']))
        print('\tDC: ' + str(x[key]['gain0']))
        for index in range(len(x[key]['poles'])):
            print('\tP' + str(index) + ":", str(x[key]['poles'][index]))
        for index in range(len(x[key]['zeros'])):
            print('\tZ' + str(index) + ":", str(x[key]['zeros'][index]))
    return None


def print_symbolic_equations(eq_list):
    """Print symbolic equations for visual inspection

    **Parameters:**

    eq_list : list
        The list of equations to be printed. This is what ``sympy`` will be
        asked to solve, typically.
    """
    print("+--")
    for eq in eq_list:
        print("| " + str(eq))
    print("+--")
    return


def print_result_check(badvars, verbose=2):
    """Prints out the results of an OP check

    It assumes one set of results is calculated with :math:`G_{min}`, the other
    without.

    **Parameters:**

    badvars : list
        The list returned by :func:`results.op_solution.gmin_check`.
    verbose : int, optional
        The verbosity level, from 0 (silent) to 6.

    """
    if len(badvars):
        print("Warning: solution is heavvily dependent on gmin.")
        print("Affected variables:")
        for bv in badvars:
            print(bv)
    else:
        if verbose:
            print("Difference check within margins.")
            print("(Voltage: er=" + str(options.ver) + ", ea=" + str(options.vea) + \
                ", Current: er=" + \
                str(options.ier) + ", ea=" + str(options.iea) + ")")
    return None


def table(data, *args, **argsd):
    """Format a fixed width table for pretty printing

    No data processing is done here, instead we call `tabulate
    <https://pypi.python.org/pypi/tabulate/>`_'s ``tabulate.tabulate()``,
    passing all arguments unmodified.

    **Parameters:**

    data : list-of-lists or a dictionary of iterables or a 2D NumPy array (or more).
        The tabular data.

    The remaining arguments, not documented here, are:

    headers : sequence, optional
        An explicit list of column headers.
    tablefmt : str, optional
        Table formatting specification.
    floatfmt : str, optional
        Floats formatting specification.
    numalign : str, optional
        Alignment flag for numbers.
    stralign : str, optional
        Alignment specification for strings, eg. "right".
    missingval : str, optional
        Element for the missing values.

    """
    return _tabulate.tabulate(data, *args, **argsd)

@contextlib.contextmanager
def printoptions(*args, **kwargs):
    """A context manager for ``numpy.set_printoptions``"""
    original = np.get_printoptions()
    np.set_printoptions(*args, **kwargs)
    yield
    np.set_printoptions(**original)

locale = os.getenv('LANG')
if not locale:
    print_warning('Locale appears not set! please export LANG="en_US.UTF-8" or'
                  ' equivalent, ')
    print_warning('or ahkab\'s unicode support is broken.')

