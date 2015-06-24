# -*- coding: iso-8859-1 -*-
# results.py
# Results module
# Copyright 2011-2013 Giuseppe Venturini

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
This module provides classes for easy, dictionary-like access to simulation
results.

Simulation results are typically returned upon successful simulation of a
circuit and the user is not expected to use their constructor, but rather
to use the methods they provide to access their data set.

Overview of the data interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The solution classes define special methods according to their simulation
type but they all subclass :class:`solution`, which provides the shared
data interface.

The interface allows for accessing the values as::


    >>> ac_sol.keys()
    ['f', 'Vn1', 'Vn2', 'I(V1)', 'I(L1)', 'I(L2)']

Where ``ac_sol`` is a generic example instance of :class:`ac_solution`.

Checking with the ``in`` construct::

    >>> 'Vn1' in ac_sol
    True

Access any variable in the solution object::

    >>> ac_sol['f']
    array([ 6098.38572827,  6102.08394991,  6105.78441425,  6109.48712265,
            6113.19207648,  6116.89927708,  6120.60872583,  6124.32042408,

            [... omissis ...]

            6463.83880528,  6467.75864729,  6471.68086639])

Iterate over the results::

    >>> for var in ac_sol:
    ...     # do something with ac_sol[var]
    ...     pass

Convenience methods are available to identify and access the independent,
swept variable, when it is available::


    >>> ac_sol.get_xlabel()
    'f'
    >>> ac_sol.get_x()
    array([ 6098.38572827,  6102.08394991,  6105.78441425,  6109.48712265,
            6113.19207648,  6116.89927708,  6120.60872583,  6124.32042408,

            [... omissis ...]

            6463.83880528,  6467.75864729,  6471.68086639])

Index of the solution classes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::

    ac_solution
    dc_solution
    op_solution
    pss_solution
    pz_solution
    symbolic_solution
    tran_solution

Module reference
~~~~~~~~~~~~~~~~

"""

from __future__ import (unicode_literals, absolute_import,
                        division, print_function)

import sys
import time
import pickle
import re

import numpy as np

from . import circuit
from . import devices
from . import printing
from . import options
from . import constants
from . import csvlib

from .py3compat import text_type
from .__version__ import __version__

csvlib.SEPARATOR = "\t"

class _mutable_data(object):
    def __init__(self):
        self._init_file_done = False
    def _add_data(self, data):
        """Add the data matrix to the results set."""
        csvlib.write_csv(self.filename, data, self.variables, append=self._init_file_done)
        self._init_file_done = True


class solution(object):
    """Base class storing a set of generic simulation results.

    This class is not meant to be accessed directly, rather it is
    subclassed by the classes for the specific simulation solutions.

    **Parameters:**

    circ : circuit instance
        the circuit instance of the simulated circuit.
    outfile : string
        the filename of the save file
    """
    def __init__(self, circ, outfile):
        self.timestamp = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
        self.netlist_file = circ.filename
        self.netlist_title = circ.title
        self.vea = options.vea
        self.ver = options.ver
        self.iea = options.iea
        self.ier = options.ier
        self.gmin = options.gmin
        self.cmin = options.cmin
        self.temp = constants.T
        self.filename = outfile
        self._init_file_done = False

        self.skip_nodes_list = []   # nodi da saltare, solo interni
        self.variables = []
        self.units = case_insensitive_dict()
        self.iter_index = 0
        # Please redefine this sol_type in the subclasses
        self.sol_type = None

    def asarray(self):
        """Return all data.

        .. note::

            This method loads to memory a possibly huge data matrix.

        """
        data, _, _, _ = csvlib.load_csv(self.filename, load_headers=[],
                                        verbose=0)
        return data

    # Access as a dictionary BY VARIABLE NAME:
    def __len__(self):
        """Get the number of variables in the results set."""
        return len(self.variables)

    def __getitem__(self, name):
        """Get a specific variable, as from a dictionary."""
        # data, headers, pos, EOF = csvlib.load_csv(...)
        try:
            data, _, _, _ = csvlib.load_csv(self.filename, load_headers=[name],
                                            nsamples=None, skip=0, verbose=0)
        except ValueError:
            raise KeyError(name)
        return data.reshape((-1,))

    def get(self, name, default=None):
        """Get a solution by variable name."""
        try:
            # data, headers, pos, EOF = csvlib.load_csv(...)
            data, _, _, _ = csvlib.load_csv(self.filename, load_headers=[name],
                                            nsamples=None, skip=0, verbose=0)
        except ValueError:
            return default
        return data.reshape((-1,))

    def has_key(self, name):
        """Determine whether the result set contains a variable."""
        return name.upper() in [v.upper() for v in self.variables]

    def __contains__(self, name):
        """Determine whether the result set contains a variable."""
        return name.upper() in [v.upper() for v in self.variables]

    def keys(self):
        """Get all of the results set's variables names."""
        return self.variables

    def values(self):
        """Get all of the results set's variables values."""
        # data, headers, pos, EOF = csvlib.load_csv(...)
        data, _, _, _ = csvlib.load_csv(self.filename,
                                        load_headers=self.variables,
                                        nsamples=None, skip=0, verbose=0)
        values = [data[i, :] for i in range(data.shape[0])]
        return values

    def items(self):
        return list(zip(self.keys(), self.values()))

    # iterator methods
    def __iter__(self):
        self.iter_index = 0
        self.iter_data, self.iter_headers, _, _ = csvlib.load_csv(self.filename)
        return self

    def next(self):
        return self.__next__()

    def __next__(self):
        if self.iter_index == len(self.iter_headers):
            self.iter_index = 0
            raise StopIteration
        else:
            next_i = self.iter_headers[self.iter_index]
            if hasattr(self.iter_data, 'shape'):
                next_d = self.iter_data[self.iter_index, :]
            else:
                next_d = self.iter_data[self.iter_index]
            nxt = next_i, next_d
            self.iter_index += 1
        return nxt

class op_solution(solution, _mutable_data):
    """OP results

    **Parameters:**

    x : ndarray
        the result set
    error : ndarray
        the residual error after solution,
    circ : circuit instance
        the circuit instance of the simulated circuit
    outfile: str
        the file to write the results to.
        Use "stdout" to write to std output.
    iterations, int, optional
        The number of iterations needed for convergence, if known.
    """
    def __init__(self, x, error, circ, outfile, iterations=0):
        solution.__init__(self, circ, outfile)
        self.sol_type = "OP"
        self.iterations = iterations

        # We have mixed current and voltage results
        # per primi vengono tanti valori di tensioni quanti sono i nodi del circuito meno
        # uno, quindi tante correnti quanti sono gli elementi definiti in tensione presenti
        # (per questo, per misurare una corrente, si può fare uso di generatori di tensione
        # da 0V)

        nv_1 = circ.get_nodes_number() - 1 # numero di soluzioni di tensione (al netto del ref)
        self.results = case_insensitive_dict()
        self.errors = case_insensitive_dict()
        self.x = x

        for index in range(nv_1):
            varname = ("V" + str(circ.nodes_dict[index + 1])).upper()
            self.variables += [varname]
            self.results.update({varname: x[index, 0]})
            self.errors.update({varname: error[index, 0]})
            self.units.update({varname: "V"})
            if circ.is_int_node_internal_only(index+1):
                self.skip_nodes_list.append(index)

        index = nv_1 - 1
        for elem in circ:
            if circuit.is_elem_voltage_defined(elem):
                index = index + 1
                varname = ("I("+elem.part_id.upper()+")").upper()
                self.variables += [varname]
                self.results.update({varname: x[index, 0]})
                self.errors.update({varname: error[index, 0]})
                self.units.update({varname: "A"})

        self._op_keys, self._op_info, self.tot_power = self._get_elements_op(circ, x)

    def __str__(self):
        str_repr = \
            (("OP simulation results for '%s'" % (self.netlist_title,)) +
            ('(netlist %s)'%(self.netlist_file,) if self.netlist_file else '') +
            ('.\nRun on %s, data file %s.\n' % \
             (self.timestamp, self.filename)))
        return str_repr + self.get_table_array()

    def __getitem__(self, name):
        """Get a specific variable, as from a dictionary."""
        if not name.upper() in [v.upper() for v in self.variables]:
            raise KeyError
        his = csvlib.get_headers_index(self.variables, [name], verbose=0)
        return self.x[his]

    def get(self, name, default=None):
        """Get a solution by variable name."""
        try:
            data = self.__getitem__(name)
        except KeyError:
            return default
        return data

    def asarray(self):
        """Get all data as a ``numpy`` array"""
        return self.x

    def get_table_array(self):
        headers = ("Variable", "Units", "Value", "Error", "%")
        table = []
        for v in self.variables:
            if self.results[v] != 0:
                relerror = self.errors[v]/self.results[v]*100.0
            else:
                relerror = 0.0
            line = (v, self.units[v], self.results[v], self.errors[v], '%d' %
                    relerror)
            table.append(line)
        return printing.table(table, headers=headers)

    def _get_elements_op(self, circ, x):
        """Returns"""
        tot_power = 0
        i_index = 0
        nv_1 = circ.get_nodes_number() - 1
        op_info = {}
        op_keys = {}
        for elem in circ:
            ports_v_v = []
            if hasattr(elem, "get_op_info"):
                if elem.is_nonlinear:
                    # build the drive ports vector
                    oports = elem.get_output_ports()
                    for index in range(len(oports)):
                        dports = elem.get_drive_ports(index)
                        ports_v = []
                        for port in dports:
                            tempv = 0
                            if port[0]:
                                tempv = x[port[0]-1]
                            if port[1]:
                                tempv = tempv - x[port[1]-1]
                            ports_v.append(tempv)
                    ports_v_v.append(ports_v)
                else:
                    port = (elem.n1, elem.n2)
                    tempv = 0
                    if port[0]:
                        tempv = x[port[0]-1]
                    if port[1]:
                        tempv = tempv - x[port[1]-1]
                    ports_v_v = ((tempv,),)
                if circuit.is_elem_voltage_defined(elem):
                    i = circ.find_vde_index(elem.part_id)
                    nv_1 = circ.get_nodes_number() - 1
                    opk, opi = elem.get_op_info(ports_v_v, x[nv_1 + i])
                else:
                    opk, opi = elem.get_op_info(ports_v_v)
                if elem.part_id[0].upper() != 'M':
                    if elem.part_id[0].upper() in op_info:
                        op_info.update({elem.part_id[0].upper():
                                        op_info[elem.part_id[0].upper()]+[opi]})
                        #assert set(opk) == set(op_keys[elem.part_id[0].upper()])
                    else:
                        op_info.update({elem.part_id[0].upper():[opi]})
                        op_keys.update({elem.part_id[0].upper():[opk]})
                else:
                    op_info.update({elem.part_id.upper():opi})
                    op_keys.update({elem.part_id.upper():[[]]})

            if isinstance(elem, devices.GISource):
                v = 0
                v = v + x[elem.n1-1] if elem.n1 != 0 else v
                v = v - x[elem.n2-1] if elem.n2 != 0 else v
                vs = 0
                vs = vs + x[elem.n1-1] if elem.sn1 != 0 else vs
                vs = vs - x[elem.n2-1] if elem.sn2 != 0 else vs
                tot_power = tot_power - v*vs*elem.alpha
            elif isinstance(elem, devices.ISource):
                v = 0
                v = v + x[elem.n1-1] if elem.n1 != 0 else v
                v = v - x[elem.n2-1] if elem.n2 != 0 else v
                tot_power = tot_power - v*elem.I()
            elif isinstance(elem, devices.VSource) or \
                 isinstance(elem, devices.EVSource):
                v = 0
                v = v + x[elem.n1-1] if elem.n1 != 0 else v
                v = v - x[elem.n2-1] if elem.n2 != 0 else v
                tot_power = tot_power - v*x[nv_1 + i_index, 0]
                i_index = i_index + 1
            elif isinstance(elem, devices.FISource):
                local_i_index = 0
                found_source = False
                for e in circ:
                    if circuit.is_elem_voltage_defined(e):
                        if isinstance(e, devices.VSource) and e.part_id.lower() == elem.source_id.lower():
                            found_source = True
                            break
                        else:
                            local_i_index += 1
                if not found_source:
                    raise RuntimeError("Sensing voltage source %s for %s not found. BUG!" %
                                       (elem.source_id, elem.part_id))
                v = 0.
                v = v + x[elem.n1 - 1] if elem.n1 != 0 else v
                v = v - x[elem.n2 - 1] if elem.n2 != 0 else v
                tot_power = tot_power - v * elem.alpha * x[nv_1 + local_i_index, 0]
            elif isinstance(elem, devices.HVSource):
                try:
                    local_i_index = circ.find_vde_index(elem.source_id)
                except ValueError:
                    raise RuntimeError("Sensing voltage source %s for %s not found. BUG!" %
                                       (elem.source_id, elem.part_id))
                local_i_index2 = circ.find_vde_index(elem.part_id)
                tot_power = tot_power - elem.alpha*x[nv_1 + local_i_index, 0]* \
                                        x[nv_1 + local_i_index2, 0]
            elif circuit.is_elem_voltage_defined(elem):
                i_index = i_index + 1

        #op_info.append("TOTAL POWER: %e W\n" % (tot_power,))
        return op_keys, op_info, tot_power

    def write_to_file(self, filename=None):
        if filename is None and self.filename is None:
            # maybe warn the user here?
            return
        if filename is None:
            filename = self.filename
        if filename != 'stdout':
            fp = printing.open_utf8(filename+"info")
        else:
            fp = sys.stdout
        fp.write(self.timestamp+"\n")
        fp.write("ahkab v. "+__version__+" (c) 2006-2015 Giuseppe Venturini\n\n")
        fp.write("Operating Point (OP) analysis\n\n")
        fp.write("Netlist: %s\nTitle: %s\n" % (self.netlist_file, self.netlist_title))
        fp.write("At %.2f K\n" % (self.temp,))
        fp.write("Options:\n\tvea = %e\n\tver = %f\n\tiea = %e\n\tier = %f\n\tgmin = %e\n" \
                 % (self.vea, self.ver, self.iea, self.ier, self.gmin))
        fp.write("\nConvergence reached in %d iterations.\n" % (self.iterations,))
        fp.write("\n========\n")
        fp.write("RESULTS:\n")
        fp.write("========\n\n")
        vtable = self.get_table_array()
        fp.write(vtable+'\n')
        fp.write("\n========================\n")
        fp.write("ELEMENTS OP INFORMATION:\n")
        fp.write("========================\n\n")
        for k in sorted(self._op_info.keys()):
            t = printing.table(self._op_info[k], headers=self._op_keys[k][0])
            fp.write(t + '\n\n')
        fp.write('Total power dissipation: %g W\n\n' % self.tot_power)
        fp.flush()
        if filename != 'stdout':
            fp.close()
            # save to .op file
            self._add_data(self.x)

    def print_short(self):
        """Print a short, essential representation of the OP results"""
        table = []
        line = []
        for v in self.variables:
            line.append("%s: %g %s" % \
                        (v, self.results[v], self.units[v]))
            if len(line) == 5:
                table.append(line)
                line = []
        if len(line) > 0: # add the last line
            line += [""]*(5 - len(line))
            table.append(line)
        print(printing.table(table))

    @staticmethod
    def gmin_check(op2, op1):
        """Checks the differences between two sets of OP results.

        It is assumed that one set of results is calculated with Gmin, the other without.

        **Parameters:**

        op1, op2: op_solution instances
            the results vectors, interchangeable

        **Returns:**

        test_fail_variables : list
            The list of the variables that did not pass the test. They are extracted from
            the op_solution objects. If the check was passed, this is an empty list.
        """

        check_failed_vars = []
        for v in op1.variables:
            abserr = abs(op2.results[v] - op1.results[v])
            if op1.units[v] == 'V':
                if abserr > options.ver*max(abs(op1.results[v]),
                                            abs(op2.results[v])) + options.vea:
                    check_failed_vars.append(v)
            elif op1.units[v] == 'A':
                if abserr > options.ier*max(abs(op1.results[v]),
                                            abs(op2.results[v])) + options.iea:
                    check_failed_vars.append(v)
            else:
                print("Unrecognized unit... Bug.")
        return check_failed_vars

    def values(self):
        """Get all of the results set's variables values."""
        return np.squeeze(self.x).tolist()

    def items(self):
        vlist = []
        for j in range(self.x.shape[0]):
            vlist.append(self.x[j, 0])
        return list(zip(self.variables, vlist))

    # iterator methods
    def __iter__(self):
        self._iter_index = 0
        return self

    def next(self):
        if self._iter_index == len(self.variables):
            self._iter_index = 0
            raise StopIteration
        else:
            nxt = self.variables[self._iter_index], \
                   self.x[self._iter_index]
            self._iter_index += 1
        return nxt

    def __next__(self):
        return self.next()


class ac_solution(solution, _mutable_data):
    """AC results

    **Parameters:**

    circ : circuit instance
        the circuit instance of the simulated circuit
    start : float
       the AC sweep frequency start value, in Hz.
    stop : float
       the AC sweep frequency stop value, in Hz.
    points : int
       the AC sweep total points.
    stype : str
       the type of sweep, ``"LOG"``, ``"LIN"`` or arb. ``"POINTS"``.
    op : op_solution
       the linearization Operating Point used to compute the results.
    outfile: str
        the file to write the results to.  Use ``"stdout"`` to write to the
        standard output.
    """
    def __init__(self, circ, start, stop, points, stype, op, outfile):
        solution.__init__(self, circ, outfile)
        self.sol_type = "AC"
        self.linearization_op = op
        self.stype = stype
        self.ostart, self.ostop, self.opoints = start, stop, points

        self.variables += ["f"]
        self.units.update({"f": "Hz"})
        self.csv_headers = [self.variables[0]]

        nv_1 = circ.get_nodes_number() - 1 # numero di soluzioni di tensione (al netto del ref)
        for index in range(nv_1):
            varname = "V%s" % str(circ.nodes_dict[index + 1])
            self.variables += [varname]
            self.units.update({varname: "V"})
            if circ.is_int_node_internal_only(index+1):
                self.skip_nodes_list.append(index)

        for elem in circ:
            if circuit.is_elem_voltage_defined(elem):
                varname = "I(%s)" % elem.part_id.upper()
                self.variables += [varname]
                self.units.update({varname: "A"})

        for i in range(1, len(self.variables)):
            self.csv_headers.append("|%s|" % self.variables[i])
            self.csv_headers.append("arg(%s)" % self.variables[i])

    def _add_data(self, data):
        """Remember to call this method with REAL data - already split in ABS and PHASE."""
        csvlib.write_csv(self.filename, data, self.csv_headers, append=self._init_file_done)
        self._init_file_done = True

    def __str__(self):
        return ("<AC simulation results for '%s' (netlist %s). %s sweep, " +
                "from %g to %g rad/sec, %d points. Run on %s, data file " +
                "%s>") % (self.netlist_title, self.netlist_file, self.stype,
                          self.ostart, self.ostop, self.opoints, self.timestamp,
                          self.filename)

    def add_line(self, frequency, x):
        frequency = np.array([[frequency]])
        xsplit = np.zeros((x.shape[0]*2, 1))
        for i in range(x.shape[0]):
            xsplit[2*i, 0] = np.abs(x[i, 0])
            xsplit[2*i+1, 0] = np.angle(x[i, 0], deg=options.ac_phase_in_deg)

        data = np.concatenate((frequency, xsplit), axis=0)
        self._add_data(data)

    def get_x(self):
        return self[self.variables[0]]

    def get_xlabel(self):
        return self.variables[0]

    def asarray(self):
        """Return all data as a (possibly huge) python matrix."""
        ## data, headers, pos, EOF = csvlib.load_csv()
        data, headers, _, _ = csvlib.load_csv(self.filename, load_headers=[],
                                              nsamples=None, skip=0, verbose=0)
        cplx_data = None
        cplx_headers = []
        re1 = '\\|(.*?)\\|'
        rg = re.compile(re1, re.IGNORECASE|re.DOTALL)

        for i in range(len(headers)):
            if headers[i].upper() == self.variables[0].upper():
                if cplx_data is None:
                    cplx_data = np.array(data[i, :].reshape((1, -1)), dtype=np.complex_)
                else:
                    cplx_data = np.vstack((cplx_data, data[i, :].reshape((1, -1))))
            else:
                m = rg.search(headers[i])
                if m: # we got a |VAR|
                    var = m.group(1)
                    cplx_headers.append(var)
                    match_phase = ('arg(%s)' % var).upper()
                    ip = [h.upper() for h in headers].index(match_phase)
                    if cplx_data is None:
                        cplx_data = np.array(data[i, :]*
                                             np.exp(1j*data[ip,:]).reshape((1, -1)),
                                             dtype=np.complex_)
                    else:
                        cplx_data = np.vstack((cplx_data,
                                               (data[i, :]*
                                                np.exp(1j*data[ip, :])).reshape((1, -1))))
        return cplx_data

    # Access as a dictionary BY VARIABLE NAME:
    def __getitem__(self, name):
        """Get a specific variable, as from a dictionary."""
        if name.upper() != 'F':
            headers = ['|%s|' % name, 'arg(%s)' % name]
        else:
            headers = [name]
        try:
            # data, headers, pos, EOF = csvlib.load_csv()
            data, headers, _, _ = csvlib.load_csv(self.filename,
                                                  load_headers=headers,
                                                  nsamples=None, skip=0,
                                                  verbose=0)
        except ValueError:
            # raise the correct exception
            raise KeyError(name)
        if len(headers) == 2:
            data = data[0, :] * np.exp(1j*data[1, :])
        else:
            data = data.reshape((-1,))
        return data

    def get(self, name, default=None):
        """Get a solution by variable name."""
        try:
            data = self.__getitem__(name)
        except KeyError:
            return default
        return data

    def values(self):
        """Get all of the results set's variables values."""
        data = self.asarray()
        values = [np.real_if_close(data[0, :])]
        for i in range(1, data.shape[0]):
            values.append(data[i, :])
        return values

    def items(self):
        values = self.values()
        return zip(self.variables, values)

    # iterator methods
    def __iter__(self):
        self.iter_index = 0
        self.iter_data = self.values()
        self.iter_headers = self.variables
        return self

class dc_solution(solution, _mutable_data):
    """DC results

       **Parameters:**

       circ : circuit instance
           the simulated circuit.
       start : float
           the DC sweep start value.
       stop : float
           the DC sweep stop value.
       sweepvar : str
           the swept variable ``part_id``.
       stype : str
           the type of sweep, ``"LOG"``, ``"LIN"`` or arb. ``"POINTS"``.
       outfile : str
           the filename of the file where the results will be written.
           Use ``"stdout"`` to write to std output.
    """
    def __init__(self, circ, start, stop, sweepvar, stype, outfile):
        solution.__init__(self, circ, outfile)
        self.sol_type = "DC"
        self.start, self.stop = start, stop
        self.stype = stype

        nv_1 = circ.get_nodes_number() - 1 # numero di soluzioni di tensione (al netto del ref)
        self.variables = [sweepvar]
        self.units = case_insensitive_dict()
        if self.variables[0][0] == 'V':
            self.units.update({self.variables[0]:'V'})
        if self.variables[0][0] == 'I':
            self.units.update({self.variables[0]:'A'})

        for index in range(nv_1):
            varname = "V%s" % (str(circ.nodes_dict[index + 1]),)
            self.variables += [varname]
            self.units.update({varname:"V"})
            if circ.is_int_node_internal_only(index+1):
                self.skip_nodes_list.append(index)

        for elem in circ:
            if circuit.is_elem_voltage_defined(elem):
                varname = "I(%s)" % (elem.part_id.upper(),)
                self.variables += [varname]
                self.units.update({varname:"A"})

    def __str__(self):
        return ("<DC simulation results for '%s' (netlist %s). %s sweep of" +
                " %s from %g to %g %s. Run on %s, data file %s>") % \
               (self.netlist_title, self.netlist_file, self.stype,
                self.variables[0].upper(), self.start, self.stop,
                self.units[self.variables[0]], self.timestamp, self.filename)

    def add_op(self, sweepvalue, op):
        """A DC sweep is made of a set of OP points.

        This method adds an OP solution and
        its corresponding sweep value to the results set.
        """
        sweepvalue = np.array([[sweepvalue]])
        x = op.asarray()
        data = np.concatenate((sweepvalue, x), axis=0)
        self._add_data(data)

    def get_x(self):
        return self.get(self.variables[0])

    def get_xlabel(self):
        return self.variables[0]

class tran_solution(solution, _mutable_data):
    """Transient results

    **Parameters:**

    circ : circuit instance
        the circuit instance of the simulated circuit.
    tstart : float
        the transient simulation start time.
    tstop : float
        the transient simulation stop time.
    op : op_solution instance
        the Operating Point (OP) used to start the transient analysis.
    method : str
        the differentiation method employed.
    outfile : str
        the filename of the save file.
        Use "stdout" to write to the standard output.
    """
    def __init__(self, circ, tstart, tstop, op, method, outfile):
        solution.__init__(self, circ, outfile)
        self.sol_type = "TRAN"
        self.start_op = op
        self.tstart, self.tstop = tstart, tstop
        self.method = method

        self._lock = False

        nv_1 = circ.get_nodes_number() - 1 # numero di soluzioni di tensione (al netto del ref)
        self.variables = ["T"]
        self.units.update({"T":"s"})

        for index in range(nv_1):
            varname = ("V%s" % (str(circ.nodes_dict[index + 1]),)).upper()
            self.variables += [varname]
            self.units.update({varname:"V"})
            if circ.is_int_node_internal_only(index+1):
                self.skip_nodes_list.append(index)

        for elem in circ:
            if circuit.is_elem_voltage_defined(elem):
                varname = ("I(%s)" % (elem.part_id.upper(),)).upper()
                self.variables += [varname]
                self.units.update({varname:"A"})

    def __str__(self):
        return ("<TRAN simulation results for '%s' (netlist %s), from %g s to" +
                " %g s. Diff. method %s. Run on %s, data file %s>") % \
               (self.netlist_title, self.netlist_file, self.tstart, self.tstop,
                self.method, self.timestamp, self.filename)

    def add_line(self, time, x):
        """This method adds a solution and its corresponding time value to the results set.
        """
        if not self._lock:
            time = np.array([[time]])
            data = np.concatenate((time, x), axis=0)
            self._add_data(data)
        else:
            raise RuntimeError("Attempting to add values to a complete " +
                               "result set.")

    def lock(self):
        self._lock = True

    def get_x(self):
        return self.get(self.variables[0])

    def get_xlabel(self):
        return self.variables[0]


class pss_solution(solution, _mutable_data):
    """PSS results

    **Parameters:**

    circ : circuit instance
        the circuit instance of the simulated circuit.
    method : str
        the PSS algorithm employed.
    period : float
        the solution period.
    outfile : str
        the filename of the save file.
        Use "stdout" to write to the std output.

    .. note::

        Instantiating ``pss_solution`` creates an *empty* data set. Call
        :func:`set_results` to initialize its data.

    """
    def __init__(self, circ, method, period, outfile):
        solution.__init__(self, circ, outfile)
        self.sol_type = "PSS"
        self.period = period
        self.method = method

        # We have mixed current and voltage results
        nv_1 = circ.get_nodes_number() - 1 # numero di soluzioni di tensione (al netto del ref)
        self.variables = ["T"]
        self.units.update({"T":"s"})

        for index in range(nv_1):
            varname = "V%s" % (str(circ.nodes_dict[index + 1]),)
            self.variables += [varname]
            self.units.update({varname:"V"})
            if circ.is_int_node_internal_only(index+1):
                self.skip_nodes_list.append(index)

        for elem in circ:
            if circuit.is_elem_voltage_defined(elem):
                varname = "I(%s)" % (elem.part_id.upper(),)
                self.variables += [varname]
                self.units.update({varname:"A"})

    def __str__(self):
        return ("<PSS simulation results for '%s' (netlist %s), period %g s. " +
                "Method: %s. Run on %s, data file %s>") % \
               (self.netlist_title, self.netlist_file, self.period, self.method,
                self.timestamp, self.filename)

    def set_results(self, t, x):
        """Set the results in the data set

        .. note::

            * All the data are set at the same time for a PSS results set.
            * Instantiating ``pss_solution`` creates an empty data set.
            * This method should be called as soon as the data is available.

        **Parameters:**

        t : ndarray
            The time. The array should be 2D with shape ``(1, N)``.
        x : ndarray
            The data corresponding to the variables.
            The array should be 2D with shape ``(M, N)``, where ``M`` is the
            number of variables in the data set.

        """
        time = np.array(t)
        data = np.concatenate((time, x), axis=0)
        self._add_data(data)

    def asarray(self):
        allvalues, _, _, _ = csvlib.load_csv(self.filename, load_headers=[],
                                             nsamples=None, skip=0, verbose=0)
        return allvalues

    def get_x(self):
        return self.get(self.variables[0])

    def get_xlabel(self):
        return self.variables[0]

class symbolic_solution(object):
    """Symbolic results

    **Parameters:**

    results_dict : dict
        the results dict returned by ``sympy.solve()``,
    substitutions : dict
        the substitutions (dictionary) employed before solving,
    circ : circuit instance
        the circuit instance of the simulated circuit.
    outfile : str, optional
        the filename of the save file.
        Use ``"stdout"`` to write to the standard output.
    tf : bool, optional
        Transfer function flag: set this to ``True`` if this set of results
        corrsponds to a transfer function.  Defaults to ``False``.
    """
    def __init__(self, results_dict, substitutions, circ, outfile=None, tf=False):
        self.sol_type = "Symbolic"
        self.timestamp = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
        self.netlist_file = circ.filename
        self.netlist_title = circ.title
        self.substitutions = substitutions
        self.tf = tf

        # the keys are strings
        # self.symbols = map(str, results_dict.keys())
        self.results = case_insensitive_dict()
        for symbol, result in results_dict.items():
            self.results.update({str(symbol):result})

        self._symbols = list(results_dict.keys()) # keep them, they're useful
        for expr in list(results_dict.values()):
            if tf:
                expr = expr['gain']
            for symb in expr.atoms():
                if symb.is_Symbol and symb not in self._symbols:
                    self._symbols.append(symb)
        self.filename = outfile if outfile != 'stdout' else None
        if self.filename is not None:
            self.save()

    def as_symbol(self, variable):
        """Converts a string to the corresponding symbolic variable.

        This symbol may then be used by the user as an atom to construct
        new expressions, modify the results expressions or it can be passed
        to Sympy's functions.

        **Parameters:**

        variable : string
            The string that identifies the variable. Eg. ``'R1'`` for the variable
            corresponding to the resistance of the resistor ``R1``. Note that the
            case is disregarded and that the first letter defines the type of
            the element (resistor, capacitor...).

        **Returns:**

        symbol : Sympy symbol
            The corresponding symbol, if it exists in the result set.

        **Raises:**

        ValueError : exception
            In case no such symbol is found.
        """

        symbs = [x for x in self._symbols if x.name.lower() == variable.lower()]
        if len(symbs) == 0:
            raise ValueError("No symbol %s in the results set."%(variable,))
        else:
            return symbs[0]

    def as_symbols(self, spacedstr):
        """Convenience function to call :func:`as_symbol` multiple times.

        **Parameters:**

        spacedstr : string,
            A string containing several symbol identifiers separated by spaces.
            Eg. ``'R1 C2 L3'``.

        **Returns:**

        (s1, s2, ...) : tuple of Sympy symbol instances
            The symbols corresponding to the identifiers in the string supplied,
            ordered as the identifiers in the string.

        **Raises:**

        ValueError : exception
            In case any corresponding symbol is not found.
        """
        return list(map(self.as_symbol, spacedstr.split()))

    def save(self):
        """Write the results to disk.

        It is necessary first to set the ``filename`` attribute, indicating
        which file to write to.

        **Raises:**

        RuntimeError : exception
            If the `filename` attribute is not set.
        """
        if not self.filename:
            raise RuntimeError("Writing the results to file requires setting the \
                              'filename' attribute")
        with open(self.filename, 'wb') as fp:
            pickle.dump(self, fp, protocol=2)

    @staticmethod
    def load(filename):
        """Static method to load a symbolic solution from disk.

        **Parameters:**

        filename : str
            The filename corresponding to the file to load from.

        **Returns:**

        sol : symbolic solution instance
            The solution instance loaded from disk.

        .. warning::

            This method employs ``pickle.load``, which is to be used exclusively
            on trusted data. **Only load trusted simulation files!**

        """
        with open(filename, 'rb') as fp:
            asolution = pickle.load(fp)
        return asolution

    def __str__(self):
        str_repr = "Symbolic %s results for '%s' (netlist %s).\nRun on %s.\n" % \
                   ('simulation'*(not self.tf) + 'transfer function'*self.tf,
                    self.netlist_title, self.netlist_file, self.timestamp)
        keys = list(self.results.keys())
        keys.sort(key=str)
        if not self.tf:
            for key in keys:
                str_repr +=  str(key) + "\t = " + str(self.results[key]) + "\n"
        else:
            for key in keys:
                str_repr +=  str(key) + ":\n\t%s:\t%s\n" % \
                             ('gain', self.results[key]['gain'])
                if self.results[key]['gain'] == self.results[key]['gain0']:
                    continue
                str_repr +=  "\t%s:\t%s\n" % \
                             ('gain0', self.results[key]['gain0'])
                for sing in ('poles', 'zeros'):
                    if not len(self.results[key][sing]):
                        continue
                    str_repr += "\t%s:\n" % sing
                    for p in self.results[key][sing]:
                        str_repr +=  "\t\t" + str(p) + "\n"
        return str_repr

    # Access as a dictionary:
    def __len__(self):
        """Get the number of variables in the results set."""
        return len(self.results)

    def __getitem__(self, name):
        """Get a specific header, as from a dictionary."""
        return self.results[str(name).upper()]

    def get(self, name, default=None):
        """Get the solution corresponding to a variable."""
        name = str(name).upper()
        try:
            return self.results[name]
        except KeyError:
            return default

    def has_key(self, name):
        """Determine whether the result set contains a variable."""
        return str(name).upper() in self.results

    def __contains__(self, name):
        """Determine whether the result set contains a variable."""
        return str(name).upper() in self.results

    def keys(self):
        """Get all of the results set's variable's names."""
        return list(self.results.keys())

    def values(self):
        """Get all of the results set's variable's values."""
        return list(self.results.values())

    def items(self):
        """Get all solutions."""
        return list(self.results.items())

    # iterator methods
    def __iter__(self):
        self.iter_index = -1
        return self

    def __next__(self):
        """Iterator method."""
        if self.iter_index == len(list(self.results.keys())) - 1:
            self.iter_index = 0
            raise StopIteration
        else:
            self.iter_index += 1
        return list(self.results.keys())[self.iter_index], \
               self.results[self._symbols[self.iter_index]]

    def next(self):
        return self.__next__()

class pz_solution(solution, _mutable_data):
    """PZ results

    **Parameters:**

    circ : circuit instance
        the circuit instance of the simulated circuit.
    poles : sequence
        the circuit zeros
    zeros : sequence
        the circuit poles
    outfile : str
        the filename of the save file.
    """
    def __init__(self, circ, poles, zeros, outfile):
        solution.__init__(self, circ, outfile)
        self.sol_type = "PZ"
        self.poles = np.sort_complex(np.array(poles).reshape((-1,)))
        self.zeros = np.sort_complex(np.array(zeros).reshape((-1,)))
        data = np.vstack((self.poles.reshape((-1, 1)),
                          self.zeros.reshape((-1, 1))))
        if np.prod(self.poles.shape):
            for i in range(self.poles.shape[0]):
                self.variables += ['p%d' % i]
        if np.prod(self.zeros.shape):
            for i in range(self.zeros.shape[0]):
                self.variables += ['z%d' % i]
            for v in self.variables:
                self.units.update({v: "rad/s"})
            self.csv_headers = []
            for i in range(len(self.variables)):
                self.csv_headers.append("Re(%s)" % self.variables[i])
                self.csv_headers.append("Im(%s)" % self.variables[i])

            # save in Re/Im form
            sdata = data.reshape(-1).view(np.float_).reshape((-1, 1))
            self._add_data(sdata)

            # store local data too:
            self.data = case_insensitive_dict()
            for i in range(len(self.variables)):
                self.data.update({self.variables[i]: data[i, 0]})

    def _add_data(self, data):
        """Remember to call this method with REAL data - already split in RE and IM."""
        csvlib.write_csv(self.filename, data, self.csv_headers, append=self._init_file_done)
        self._init_file_done = True

    def __str__(self):
        return ("PZ simulation results for %s (netlist %s).\n" + \
               "Poles: %s\nZeros: %s") % \
               (self.netlist_title, self.netlist_file, list(self.poles),
                list(self.zeros))

    # Access as a dictionary BY VARIABLE NAME:
    def __getitem__(self, name):
        """Get a specific variable, as from a dictionary."""
        if name in self.data:
            return self.data[name]
        raise KeyError

    def get(self, name, default=None):
        try:
            return self.data[name]
        except KeyError:
            return default

    def has_key(self, name):
        """Determine whether the result set contains a variable."""
        return name in self.data

    def __contains__(self, name):
        """Determine whether the result set contains a variable."""
        return name in self.data

    def keys(self):
        """Get all of the results set's variable's names."""
        return self.data.keys()

    def values(self):
        """Get all of the results set's variable's values."""
        return self.data.values()

    def items(self):
        return self.data.items()

    # iterator methods
    def __iter__(self):
        # take into account that we increment first, then return
        # the value
        self.iter_index = -1
        return self

    def __next__(self):
        # redefine this or the superclass method will be used on
        # PY3, calling the superclass' next() and not *our* next
        # method below.
        return self.next()

    def next(self):
        if self.iter_index == len(self.variables)-1:
            self.iter_index = 0
            raise StopIteration
        else:
            self.iter_index += 1
        return self.variables[self.iter_index], \
               self.data[self.variables[self.iter_index]]

class case_insensitive_dict(object):
    """A dictionary that uses case-insensitive strings as keys.
    """
    def __init__(self):
        self._dict = {}

    def __len__(self):
        """Get the number of elements in the set."""
        return len(self._dict)

    def __str__(self):
        rpr = "{"
        keys = list(self._dict.keys())
        for i in range(len(keys)):
            rpr += "%s: %s" % (keys[i], self._dict[keys[i]])
            if i < len(keys) - 1:
                rpr += ", "
            else:
                rpr += "}"
        return rpr

    def __getitem__(self, name):
        """Get a specific variable, as from a dictionary."""
        keys = list(self._dict.keys())
        try:
            i = [k.upper() for k in keys].index(text_type(name).upper())
        except ValueError:
            raise KeyError(name)
        return self._dict[keys[i]]

    def get(self, name, default=None):
        """Given the case-insensitive string key ``name``, return its corresponding value.

        If not found, return ``default``.
        """
        try:
            keys = list(self._dict.keys())
            i = [k.upper() for k in keys].index(name.upper())
        except ValueError:
            return default
        return self._dict[keys[i]]

    def has_key(self, name):
        """Determine whether the result set contains the variable ``name``."""
        return name.upper() in [k.upper() for k in list(self._dict.keys())]

    def __contains__(self, name):
        """Determine whether the result set contains a variable."""
        return name.upper() in [k.upper() for k in list(self._dict.keys())]

    def keys(self):
        """Get all keys"""
        return list(self._dict.keys())

    def values(self):
        """Get all values"""
        return list(self._dict.values())

    def items(self):
        """Get all keys and values pairs"""
        return list(self._dict.items())

    def update(self, adict):
        """Update the dictionary contents with the mapping in the dictionary ``adict``.
        """
        return self._dict.update(adict)

    # iterator methods
    def __iter__(self):
        """Iterator"""
        return self._dict.__iter__()
