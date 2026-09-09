# -*- coding: iso-8859-1 -*-
# plotting.py
# Module to plot the simulation data (through matplotlib)
# Copyright 2009-2013 Giuseppe Venturini

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
This module offers the functions needed to plot the results
of a simulation.

It is only functional if `matplotlib <http://matplotlib.org/>`_
is installed.

Module reference
''''''''''''''''

"""

from __future__ import (unicode_literals, absolute_import,
                        division, print_function)

import re
import numpy as np

from . import printing
from . import options
from . import py3compat

try:
    import pylab
except ImportError:
    # no matplotlib on a platform that supports it
    if not py3compat.PYPY:
        printing.print_warning('Matplotlib not found.')

def _split_netlist_label(labels_string):
    labels_string = labels_string.strip().upper()
    ret_labels = []
    p = re.compile(r'V\s*\(\s*(\w*)\s*,\s*(\w*)\s*\)', re.IGNORECASE)
    labels_list = p.findall(labels_string)
    for i in range(len(labels_list)):
        l2 = "V" + labels_list[i][0]
        l1 = "V" + labels_list[i][1]
        if l1 != 'V0':
            ret_labels.append((l2, l1))
        else:
            ret_labels.append((l2,))
    p = re.compile(r'V\s*\(\s*(\w*)\s*\)', re.IGNORECASE)
    labels_list = p.findall(labels_string)
    for i in range(len(labels_list)):
        l2 = "V" + labels_list[i]
        l1 = None
        ret_labels.append((l2, l1))
    p = re.compile(r'I\s*\(\s*(\w*)\s*\)', re.IGNORECASE)
    labels_list = p.findall(labels_string)
    for i in range(len(labels_list)):
        l2 = "I(" + labels_list[i] + ")"
        l1 = None
        ret_labels.append((l2, l1))
    if len(ret_labels) == 0:
        raise ValueError("Unrecognized plot labels: " + labels_string)
    return ret_labels


def _setup_plot(fig, title, xvu, yvu, log=False, xlog=False, ylog=False):
    """Setup the figure for plotting. 

    **Parameters:**

    fig : figure
        A matplotlib figure
    title : string
        The plot title:
    xvu : tuple
        A tuple defined as ``xvu = (xvar, xunit)``, where ``xvar`` is the
        x-axis variable label (str) and ``xunit`` is its unit (str too).
    yvu : list of tuples
        defined as yvu += [
        Each tuple defined as ``(yvar, yunit)``, where ``yvarN`` is the
        y-axis variable label (str) and ``yunit`` is its unit (str too).
    log : bool, optional
        Whether to set all scales to ``log``.
    xlog : bool, optional
        Whether to set the x-axis scale to ``log``.
    ylog : bool, optional
        Whether to set the y-axis scale to ``log``.

    """
    # set the currently active figure to fig
    pylab.figure(fig.number)
    # xvar, xunit = xvu
    pylab.title(title.upper())
    ax = pylab.gca()

    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    ax.xaxis.grid(False)
    ax.yaxis.grid(False)

    if log or xlog:
        ax.set_xscale('log')
    pylab.xlabel('%s [%s]' % xvu)
    yunits = []
    yinitials = []
    for yv, yu in yvu:
        yv = yv[:].replace('|', "")
        if yu not in yunits:
            yunits.append(yu)
            yinitials.append(yv[0])
    ylabel = ""
    for yi, yu in zip(yinitials, yunits):
        ylabel += "%s [%s] , " % (yi, yu)
    ylabel = ylabel[:-3]
    pylab.ylabel(ylabel)

    if log or ylog:
        ax.set_yscale('log')
    # fig.tight_layout()


def save_figure(filename, fig=None):
    """Apply the figure options for saving and then save the supplied
    figure to ``filename``.
    
    The format of the output figure is set by ``options.plotting_outtype``.

    **Parameters:**

    filename : string
        The output filename.
    fig : figure object, optional
        The figure to be saved.

    **Returns:**
    
    ``None``.
    """
    if fig is None:
        fig = pylab.gcf()
    fig.set_size_inches(*options.plotting_save_figsize)
    pylab.savefig(filename, dpi=100, bbox_inches='tight',
                  format=options.plotting_outtype)
    fig.set_size_inches(*options.plotting_display_figsize)

def _data_abs_arg_pass(res, label):
    # extract abs / phase if needed or pass the data
    if label[0] == label[-1] == '|':
        data = np.absolute(res[label[1:-1]])
        units = res.units[label[1:-1]]
    elif label[0:4] == 'arg(' and label[-1] == ')':
        data = np.angle(res[label[4:-1]], deg=options.ac_phase_in_deg)
        units = res.units[label[4:-1]]
    else:
        data = res[label]
        units = res.units[label]
    return data, units

def plot_results(title, y2y1_list, results, outfilename=None):
    """Plot the results.

    **Parameters:**

    title : string
        The plot title
    y2y1_list : list
        A list of tuples. Each tuple has to be in the format ``(y2, y1)``.
        Each member of the tuple has to be a valid identifier. You can 
        check the possible voltage and current identifiers in the 
        result set calling ``res.keys()``, where ``res`` is a solution
        object.
    result : solution object or derivate
        The results to be plotted.
    outfilename : string, optional
        The filename of the output file. If left unset, the plot will
        not be written to disk. The format is set through
        ``options.plotting_outtype``.

    **Returns:**

    ``None``.
    """
    if results is None:
        printing.print_warning("No results available for plotting. Skipping.")
        return
    fig = pylab.figure(figsize=options.plotting_display_figsize)
    analysis = results.sol_type.upper()
    gdata = []
    x, xlabel = results.get_x(), results.get_xlabel()
    xunit = results.units[xlabel]
    yvu = []

    for y2label, y1label in y2y1_list:
        if y1label is not None and y1label != '':
            try:
                data1, _ = _data_abs_arg_pass(results, y1label)
            except ValueError as e:
                printing.print_warning(str(e) + " " + y1label)
                continue
            line_label = y2label + "-" + y1label
        else:
            line_label = y2label
            data1 = 0
        try:
            data2, units = _data_abs_arg_pass(results, y2label)
        except ValueError as e:
            printing.print_warning(str(e) + " " + y2label)
            continue
        yvu += [(line_label, units)]
        gdata.append((data2 - data1, line_label))

    if xlabel == 'f':
        xlog = True
    else:
        xlog = False
    _setup_plot(fig, title, (xlabel, xunit), yvu, xlog=xlog)
    ms = 7./(1. + max(np.log(len(x)/100.), 0))
    ms = ms if ms > 2 else 0.

    ymax, ymin = None, None
    for y, label in gdata:
        [line] = pylab.plot(
            x, y, options.plotting_style, label=label +
            " (" + analysis + ")", ms=ms,
            mfc='w', lw=options.plotting_lw, mew=options.plotting_lw)
        line.set_mec(line.get_color())  # nice empty circles
        ymax = y.max() if ymax is None or y.max() > ymax else ymax
        ymin = y.min() if ymin is None or y.min() < ymin else ymin
    pylab.xlim((x.min(), x.max()))
    pylab.ylim((ymin - (ymax - ymin) * .01, ymax + (ymax - ymin) * .01))
    pylab.grid(True)
    pylab.legend()

    if outfilename is not None and options.plotting_outtype is not None:
        save_figure(outfilename, fig)
    return


def show_plots():
    """See the fruit of your work!"""
    pylab.show()
