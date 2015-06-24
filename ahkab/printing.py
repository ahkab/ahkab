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
    print_fourier
    print_spicefft
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
    if py3compat.PY2:
        fp = open(filename, 'w')
        UTF8Writer = codecs.getwriter('utf8')
        fp = UTF8Writer(fp)
    else:
        fp = open(filename, 'w', encoding='UTF-8')
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

def print_fourier(label, f, F, THD, outfile='stdout'):
    """Print the results of a Fourier postprocess

    **Parameters:**

    label : str, or tuple of str
        The identifier of a variable. Eg. ``'Vn1'`` or ``'I(VS)'``.
        If a tuple of two identifiers is provided, it will be interpreted as the
        difference of the two, in the form ``label[0]-label[1]``.
    f : ndarray of floats
        The frequencies, including the DC.
    F : ndarray of complex data
        The result of the Fourier transform.
    THD : float, optional
        The total harmonic distortion, if ``freq`` was specified in the FFT
        analysis. Defaults to ``None``.
    outfile : str, optional
        The file name to print to. Defaults to ``'stdout'``, for the standard
        output.
     """
    if type(label) in py3compat.string_types:
        pass
    elif type(label) == tuple:
        if len(label) == 1 or (len(label) >= 2 and label[1] is None):
            label = label[0]
        else:
            label = '%s - %s' % (label[0], label[1])
    if outfile == 'stdout':
        fp = sys.stdout
    else:
        fp = open(outfile, 'w')
    fp.write('Fourier components of transient response %s\n' % label)
    DC = np.real_if_close(F[0])
    fp.write('DC component =  %g\n' % DC)
    headers = ('#', 'Frequency [Hz]', 'Magnitude', 'Normalized', 'Phase [deg]',
               'Normalized [deg]')
    ffttable = []
    i = np.argmax(np.abs(F))
    ref_f = f[i]
    ref_amplitude = np.abs(F)[i]
    ref_angle = np.angle(F, deg=True)[i]
    fp.write('Normalization to f = %g Hz, |A| = %g, phase(A) = %g deg\n' %
             (ref_f, ref_amplitude, ref_angle))
    for i in range(1, len(F)):
        ffttable += [[i, f[i], np.abs(F[i]), abs(F[i])/ref_amplitude,
                      np.angle(F[i], deg=True),
                      np.angle(F[i], deg=True) - ref_angle]]
    fp.write(table(ffttable, headers=headers))
    fp.write('\nTotal harmonic distortion: %g %%\n' % (100*THD))
    if outfile != 'stdout':
        fp.close()

def print_spicefft(label, f, F, THD=None, uformat='NORM', window=None,
                   outfile='stdout'):
    """Print the results of an FFT postprocess

    **Parameters:**

    label : str, or tuple of string
        The identifier of a variable. Eg. ``'Vn1'`` or ``'I(VS)'``.
        If a tuple of two identifiers is provided, it will be interpreted as the
        difference of the two, in the form ``label[0]-label[1]``.
    f : ndarray of floats
        The frequencies, including the DC.
    F : ndarray of complex data
        The result of the Fourier transform.
    THD : float, optional
        The total harmonic distortion, if ``freq`` was specified in the FFT
        analysis. Defaults to ``None``.
    uformat : str, optional
        The parameter format selects whether normalized or unnormalized
        magnitudes are printed. It is to be set to 'NORM' (default value) for
        normalized magnitude, to 'UNORM' for unnormalized.
    window : str, optional
        The window employed in the FFT analisys. Defaults to rectangular.
    outfile : str
        The file name to print to. Defaults to ``'stdout'``, for the standard
        output.
    """
    if type(label) in py3compat.string_types:
        pass
    elif type(label) == tuple:
        if len(label) == 1 or (len(label) >= 2 and label[1] is None):
            label = label[0]
        else:
            label = '%s - %s' % (label[0], label[1])
    if uformat.upper() not in ('NORM', 'UNORM'):
        raise ValueError('fft(): format may be "NORM" or "UNORM", got %s' %
                         uformat)
    angle = np.angle(F, deg=True)
    F = 20*np.log10(abs(F))  # dB
    if uformat.upper() == 'NORM':
        i = np.argmax(F)
        fnorm = f[i]
        Fnorm = F[i]
        F -= F.max()
        unit = 'dBc'
    else:
        unit = 'dB'
    if window is not None:
        window = window.upper()
        if window not in (options.RECT_WINDOW, options.BART_WINDOW,
                               options.HANN_WINDOW, options.HAMM_WINDOW,
                               options.BLACK_WINDOW, options.HARRIS_WINDOW,
                               options.GAUSS_WINDOW, options.KAISER_WINDOW):
            raise ValueError(('fft(): window may be %s, %s, %s, %s, %s, %s, %s'+
                              ' or %s, got %s') % (options.RECT_WINDOW,
                                                   options.BART_WINDOW,
                                                   options.HANN_WINDOW,
                                                   options.HAMM_WINDOW,
                                                   options.BLACK_WINDOW,
                                                   options.HARRIS_WINDOW,
                                                   options.GAUSS_WINDOW,
                                                   options.KAISER_WINDOW,
                                                   window))
    else:
        window = options.RECT_WINDOW
    # set up the output file
    if outfile == 'stdout':
        fp = sys.stdout
    else:
        fp = open(outfile, 'w')
    fp.write('FFT components of transient response %s\n' % label)
    fp.write('Window: %s\n' % options.WINDOWS_NAMES[window.upper()])
    fp.write('Start Frequency: %g Hz\n' % f[1])
    fp.write('Stop Frequency: %g Hz\n' % f[-1])
    if uformat.upper() == 'NORM':
        fp.write(('Normalization to carrier at %g Hz: '% fnorm) +
                 ('magnitude %g dB\n' % Fnorm))
    fp.write(('DC component: mag = %g ' + unit + '\n') % F[0])
    headers = ('#', 'Frequency [Hz]', 'Magnitude ['+unit+']', 'Phase [deg]')
    ffttable = []
    for i in range(1, len(F)):
        ffttable += [[i, f[i], F[i], angle[i]]]
    fp.write(table(ffttable, headers=headers))
    fp.write('\n')
    if THD is not None:
        fp.write('Total harmonic distortion: %g %%\n' % (100*THD))
    if outfile != 'stdout':
        fp.close()

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

