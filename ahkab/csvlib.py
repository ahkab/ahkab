""" This file contains common routines for handling
Comma Separated Values (CSV) files.

Functions:

1. CSV write/load:
    write_csv(filename, data, headers, append=False)
    load_csv(filename, load_headers=[], nsamples=None, skip=0L, verbose=3)

2. MISC utilities
    get_headers_index(headers, load_headers, verbose=3):
    get_csv_headers(filename):

3. Internal routines
    _get_fp(filename, mode='r')
    _close_fp(fp, filename)

"""

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

    mode = 'a' if append else 'w'
    fp = _get_fp(filename, mode)

    if not data.shape[0] == len(headers):
        print "(W): write_csv(): data and headers don't match. Continuing anyway."
        print "DATA: " + str(data.shape) + " headers length: " + str(len(headers))

    headers = SEPARATOR.join(headers) if not append else ""
    numpy.savetxt(fp, data.T, delimiter=SEPARATOR, header=headers, comments='#')

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
    if filename == 'stdout' or filename == '-' or filename == sys.stdout:
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


def get_headers_index(headers, load_headers=[], verbose=3):
    """Creates a list of integers. Each element in the list is the COLUMN index
    of the signal according to the supplied headers.

    headers: list of strings, the signal names, as returned by get_csv_headers()

    load_headers : list, optional
        The headers for the data to be loaded. If not provided, all indeces will
        be returned.

    Returns a list of int."""
    if not len(load_headers):
        return list(range(len(headers)))
    his = []
    lowcase_headers = map(str.lower, headers)

    for lh in load_headers:
        try:
            his = his + [lowcase_headers.index(lh.lower())]
        except ValueError:
            if verbose:
                print "(W): header " + lh + " not found. Skipping."
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

    headers = get_csv_headers(filename)
    his = get_headers_index(headers, load_headers, verbose=verbose)
    if len(load_headers) and len(his) != len(load_headers):
        raise ValueError("Specified header not found")

    fp = _get_fp(filename, mode="r")
    data = numpy.loadtxt(fp, delimiter=SEPARATOR, usecols=his, unpack=True, skiprows=skip, ndmin=2)
    _close_fp(fp, filename)

    # prepare return values
    EOF = (nsamples is None) or (nsamples == data.shape[1])
    if nsamples is not None:
        data = data[:, :min(nsamples, data.shape[1])]
    pos = skip + data.shape[1]
    headers = map(headers.__getitem__, his)

    return data, headers, pos, EOF
