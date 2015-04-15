# -*- coding: utf-8 -*-
# csvlib.py
# Implementation of routines for CSV data I/O
# Copyright 2012 Giuseppe Venturini
#
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
"""The ``csvlib`` module contains common routines for handling
Comma Separated Values (CSV) or Tab Separated Values (TSV) files.

Functions:

1. CSV write/load:

  * :func:`write_csv`
  * :func:`load_csv`

2. MISC utilities

  * :func:`get_headers`
  * :func:`write_headers`
  * :func:`get_headers_index`

The separator can be selected setting:

::

   csvlib.SEPARATOR = '\\t' # default value


"""

# Additionally the following internal functions are
# available:
#3. Internal routines
#    _get_fp(filename, mode='rb')
#    _close_fp(fp, filename)

from __future__ import (unicode_literals, absolute_import,
                        division, print_function)

import io
import sys
import copy
import os
import numpy as np

from . import options

SEPARATOR = u"\t"


def write_csv(filename, data, headers, append=False):
    """Writes data in CVS format to filename.

    The headers have to be ordered according to the data order.

    **Parameters:**

    filename : string
        the path to the file to be written.
        Use 'stdout' to write to stdout

    data : ndarray
        The data to be written. Notice that variables are swept across *rows*,
        time samples are swept along *columns*.
        Or equivalently: ``data[variable_index, sample_number]``

    headers : list of strings
        the signal names, ordered so that ``headers[i]`` corresponds to 
        ``data[i, :]``.

    append : bool, optional
        If False, the file (if it exists) will be rewritten, otherwise
        it will be appended to.

    """

    mode = 'ab' if append else 'wb'
    fp = _get_fp(filename, mode)

    if not data.shape[0] == len(headers):
        print("(W): write_csv(): data and headers don't match. Continuing anyway.")
        print("DATA: " + str(data.shape) + " headers length: " + str(len(headers)))

    headers = SEPARATOR.join(headers) if not append else ""
    np.savetxt(fp, data.T, delimiter=SEPARATOR, header=headers, comments='#')

    _close_fp(fp, filename)


def write_headers(filename, headers):
    """Writes headers in CVS format to filename.

    **Parameters:**

    filename : string
        the path to the file to be written.
        Use 'stdout' to write to stdout.

    headers : a list of strings
        the signal names, ordered.

    """
    fp = _get_fp(filename, mode="wb")
    headers = copy.copy(headers)
    if not headers[0][0] == '#':
        headers[0] = '#' + headers[0]
    for hi in range(len(headers)):
        fp.write(headers[hi])
        if hi < len(headers) - 1:
            fp.write(SEPARATOR)
        else:
            fp.write("\n")


def _get_fp(filename, mode="r"):
    if filename == 'stdout' or filename == '-' or filename == sys.stdout:
        if mode == 'w' or mode == 'a' or mode == 'wb' or mode == 'ab':
            fp = sys.stdout if not hasattr(sys.stdout, 'buffer') else sys.stdout.buffer
        else:
            print("(EE) Mode %s is not supported for stdout." % (mode,))
            fp = None
    else:
        fp = io.open(filename, mode)#, encoding=options.encoding)
    return fp


def _close_fp(fp, filename):
    try:
        fp.flush()
    except IOError:
        pass
    if filename == 'stdout':
        pass
    else:
        fp.close()


def get_headers_index(headers, load_headers=None, verbose=3):
    """Creates a list of integers. Each element in the list is the COLUMN index
    of the signal according to the supplied headers.

    **Parameters:**

    headers : list of strings,
        the signal names, as returned by :func:`get_headers`.

    load_headers : list, optional
        The headers for the data to be loaded. If not provided, all indeces will
        be returned.

    **Returns:**

    the header indeces : a list of int.

    """
    if load_headers is None or not len(load_headers):
        return list(range(len(headers)))
    his = []
    lowcase_headers = [i.lower() for i in headers]

    for lh in load_headers:
        try:
            his = his + [lowcase_headers.index(lh.lower())]
        except ValueError:
            if verbose:
                print("(W): header " + lh + " not found. Skipping.")
    return his


def get_headers(filename):
    """Reads the signals inside a file.

    The order of the signals in the list corresponds to the order of the
    signals in the file.

    *Parameters:*

    filename : string
        the path to the file from which the header is to be read

    **Returns:**

    headers : list of strings.

    """

    fp = _get_fp(filename, 'r')
    headers = None
    line = ""
    while line == "":
        line = fp.readline()
        line = line.strip()
        if line[0] == '#':
            line = line[1:]
    headers = line.split(SEPARATOR)
    return headers


def load_csv(filename, load_headers=None, nsamples=None, skip=0, verbose=3):
    """Reads data in CVS format from filename.

    Supports:

    * selective signal loading,
    * loading up to a certain number of samples,
    * skipping to a certain line, to allow incremental reading of big files.

    **Parameters:**

    filename : string
        the path to the file to be read.

    load_headers : list of strings, optional
        Each one being a signal to be loaded. An empty list (or None) is 
        interpreted as "read all signals".

    nsamples : int, optional
        The number of samples to be read for each signal. If ``None``,
        read all available samples.

    skip : int, optional
        The index of the first sample to be read. Default: 0

    **Returns:**

    data : ndarray 
        The data, ordered according to the order of ``load_headers``
        (or the order on file if ``load_headers`` was empty),

    headers : list of strings
        the names of the signals read from file,

    pos : int
        position of the last sample read +1, referred to the
        sample #0 in the file.

    EOF : bool
        A flag set to true is all the data in the file were read.

    """

    if filename == 'stdout':
        print("Can't load data from stdout.")
        return None, None, None, None

    headers = get_headers(filename)
    his = get_headers_index(headers, load_headers, verbose=verbose)
    if load_headers and len(his) != len(load_headers):
        raise ValueError("Specified header not found")

    fp = _get_fp(filename)
    data = np.loadtxt(fp, delimiter=SEPARATOR, usecols=his, unpack=True, skiprows=skip, ndmin=2)
    _close_fp(fp, filename)

    # prepare return values
    EOF = (nsamples is None) or (nsamples == data.shape[1])
    if nsamples is not None:
        data = data[:, :min(nsamples, data.shape[1])]
    pos = skip + data.shape[1]
    headers = list(map(headers.__getitem__, his))

    return data, headers, pos, EOF
