""" This file contains common routines for handling
Comma Separated Values (CSV) files.

Functions:

1. CSV write/load:
    write_csv(filename, data, headers, append=False)
    load_csv(filename, load_headers=[], nsamples=None, skip=0L, verbose=3)

2. MISC utilities
    get_headers_index(headers, load_headers. verbose=3):
    get_csv_headers(filename):

3. Internal routines
    _get_fp(filename, mode='r')
    _close_fp(fp, filename)

"""

__version__ = "0.091"

import sys
import copy
import numpy

SEPARATOR = "\t"


def write_csv(filename, data, headers, append=False):
    """Writes data in CVS format to filename.
    The headers have to be ordered according to the data order.

    filename: a string, the path to the file to be written.
              use 'stdout' to write to stdout
    data: a numpy.array instance,
          variables are swept across ROWS
              time samples are swept along COLUMNS
              equivalently: data[variable_index, sample_number]
    headers: a list of strings, the signal names
             headers[i] corresponds to data[i, :]
    append: boolean value. If False, the file (if it exists)
            will be rewritten, otherwise it will be appended to.

    Returns: None
    """

    #   Writing data in CSV format to `filename`
    if not append:
        fp = _get_fp(filename, mode="w")
        headers = copy.copy(headers)
        if not headers[0][0] == '#':
            headers[0] = '#' + headers[0]
        for hi in range(len(headers)):
            fp.write(headers[hi])
            if hi < len(headers) - 1:
                fp.write(SEPARATOR)
            else:
                fp.write("\n")
    #   Appending data in CSV format to `filename`, no headers
    else:
        fp = _get_fp(filename, mode="a")

    if not data.shape[0] == len(headers):
        print "(W): write_csv(): data and headers don't match. Continuing anyway."
        print "DATA: " + str(data.shape) + " headers length: " + str(len(headers))

    for j in range(data.shape[1]):
        for i in range(len(headers)):
            fp.write("{0:g}".format(data[i, j]))
            if i < len(headers) - 1:
                fp.write(SEPARATOR)
        fp.write("\n")
    _close_fp(fp, filename)


def write_headers(filename, headers):
    """Writes headers in CVS format to filename."""
    fp = _get_fp(filename, mode="w")
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
    if filename == 'stdout':
        if mode == 'w' or mode == 'a':
            fp = sys.stdout
        else:
            print "(EE) Mode %s is not supported for stdout." % (mode,)
            fp = None
    else:
        fp = open(filename, mode)
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


def get_headers_index(headers, load_headers, verbose=3):
    """Creates a list of integers. Each element in the list is the COLUMN index
    of the signal according to the supplied headers.

    headers: list of strings, the signal names, as returned by get_csv_headers()

    Returns a list of int."""
    his = []
    lowcase_headers = map(str.lower, headers)

    for lh in load_headers:
        try:
            his = his + [lowcase_headers.index(lh.lower())]
        except ValueError:
            if verbose:
                print "(W): header " + lh + " not found. Skipping."
        except AttributeError:
            if verbose:
                print "(W): got spurious header " + str(lh) + " (not a string). Skipping."
    return his


def get_csv_headers(filename):
    """Reads the file (filename) and returns a list of the signals inside.
    The order of the signals in the list corresponds to the order of the signals in the file.

    filname: the path to the file to be open

    Returns a list of strings.
"""

    fp = _get_fp(filename, mode="r")
    headers = None
    line = ""
    while line == "":
        line = fp.readline()
        line = line.strip()
        if line[0] == '#':
            line = line[1:]
    headers = line.split(SEPARATOR)
    return headers


def load_csv(filename, load_headers=[], nsamples=None, skip=0L, verbose=3):
    """Reads data in CVS format from filename.

    Supports:
    - selective signal loading,
    - loading up to a certain number of samples,
    - skipping to a certain line.

    to allow incremental reading of big files.

    filename: a string, the path to the file to be read
    load_headers: a list of strings, each one being a signal
                  to be loaded. Empty string -> read all signals
    nsamples: int/long, number of samples to be read for each
              signal. If None, read all available samples.
    skip: index of the first sample to be read. Default: 0

    Returns:
    data: numpy.array containing the data, ordered according to
    the order of load_headers (or the order on file),
    headers: the names of the signals read from file,
    pos: long, position of the last sample read +1, referred to the
    sample #0 in the file.
    EOF: boolean, True if the EOF was reached."""

    if filename == 'stdout':
        print "Can't load data from stdout."
        return None, None, None, None

    fp = _get_fp(filename, mode="r")
    headers = None
    data = None
    sample_index = 0L
    line_index = 0
    EOF = False

    for line in fp:  # .readlines():
        line = line.strip()
        if line == '':
            continue
        if line[0] == '#' and headers is None:
            line = line[1:]
        if line[0] == '#' and headers is not None:
            continue  # comment
        if headers is None:
            headers = line.split(SEPARATOR)
            if len(load_headers):
                his = get_headers_index(headers, load_headers, verbose=verbose)
            else:
                his = range(len(headers))
        else:
            line_index = line_index + 1
            if line_index < skip:
                continue
            if data is None:
                data = numpy.zeros((len(his), 1))
            else:
                data = numpy.concatenate(
                    (data, numpy.zeros((len(his), 1))), axis=1)
            data_values = line.split(SEPARATOR)
            for i in range(len(data_values)):
                if his.count(i) > 0:
                    data[his.index(i), -1] = float(data_values[i])
                else:
                    pass
            sample_index = sample_index + 1
            if nsamples is not None and sample_index == nsamples:
                break
    _close_fp(fp, filename)

    if data is None or line == '' or nsamples > sample_index:
        EOF = True
    else:
        EOF = False

    pos = skip + sample_index

    headers = map(headers.__getitem__, his)

    return data, headers, pos, EOF
