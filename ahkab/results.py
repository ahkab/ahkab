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
This is the results module of the simulator.
"""

__version__ = "0.091"

import sys, time, copy, pickle
import numpy
import circuit, devices, printing, options, constants, csvlib
VERSION = "0.07"
csvlib.SEPARATOR = "\t"

class solution:
    def __init__(self, circ, outfile):
        """Base class storing a set of generic simulation results.
            circ: the circuit instance of the simulated circuit
            outfile: the filename of the save file
        """
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

        nv_1 = len(circ.nodes_dict) - 1 # numero di soluzioni di tensione (al netto del ref)
        self.skip_nodes_list = []   # nodi da saltare, solo interni
        self.variables = []
        self.units = case_insensitive_dict()

    def _add_data(self, data):
        """Add the data matrix to the results set."""
        csvlib.write_csv(self.filename, data, self.variables, append=self._init_file_done)
        self._init_file_done = True

    def get_type(self):
        """Please redefine this function in the subclasses."""
        raise Exception, "Undefined"

    def asmatrix(self, verbose=3):
        """Return all data as a (possibly huge) python matrix."""
        data, headers, pos, EOF = csvlib.load_csv(self.filename, load_headers=[], 
                                                  nsamples=None, skip=0L, verbose=verbose)
        return data.T

    # Access as a dictionary BY VARIABLE NAME:
    def __len__(self):
        """Get the number of variables in the results set."""
        return len(self.variables)

    def __getitem__(self, name):
        """Get a specific variable, as from a dictionary."""
        data, headers, pos, EOF = csvlib.load_csv(self.filename, load_headers=[name], 
                                                  nsamples=None, skip=0L, verbose=0)
        return data.T

    def get(self, name, default=None, verbose=3):
        """Get a solution by variable name."""
        try:
            data, headers, pos, EOF = csvlib.load_csv(self.filename, 
                                                      load_headers=[name], 
                                                      nsamples=None, skip=0L, verbose=verbose)
        except KeyError:
            return default
        return data.T

    def has_key(self, name):
        """Determine whether the result set contains a variable."""
        return name.upper() in map(str.upper, self.variables)

    def __contains__(self, name):
        """Determine whether the result set contains a variable."""
        return name.upper() in map(str.upper, self.variables)

    def keys(self):
        """Get all of the results set's variables names."""
        return self.variables

    def values(self, verbose=3):
        """Get all of the results set's variables values."""
        data, headers, pos, EOF = csvlib.load_csv(self.filename, 
                                                  load_headers=self.variables, 
                                                  nsamples=None, skip=0L, verbose=verbose)
        return data.T

    def items(self, verbose=3):
        data, headers, pos, EOF = csvlib.load_csv(self.filename, 
                                                  load_headers=self.variables, 
                                                  nsamples=None, skip=0L, verbose=verbose)
        vlist = []
        for j in range(data.shape[0]):
            vlist.append(data[j,:].T)
        return zip(headers, vlist)

    # iterator methods
    def __iter__(self):
        self.iter_index = 0
        self.iter_data, self.iter_headers, pos, EOF = csvlib.load_csv(self.filename, 
                                                                      load_headers=[], 
                                                                      nsamples=None, 
                                                                      skip=0L)
        return self

    def next(self):
        if self.iter_index == len(self.iter_headers):
            self.iter_index = 0
            raise StopIteration
        else:
            next = self.iter_headers[self.iter_index], \
                   self.iter_data[self.iter_index,:].T
            self.iter_index += 1
        return next

class op_solution(solution):
    def __init__(self, x, error, circ, outfile, iterations=0):
        """Holds a set of Operating Point results.
            x: the result set
            error: the residual error after solution,
            circ: the circuit instance of the simulated circuit
            print_int_nodes: a boolean to be set True if you wish to see voltage values 
            of the internal nodes added automatically by the simulator.
        """
        solution.__init__(self, circ, outfile)
        self.iterations = iterations

        # We have mixed current and voltage results
        # per primi vengono tanti valori di tensioni quanti sono i nodi del circuito meno 
        # uno, quindi tante correnti quanti sono gli elementi definiti in tensione presenti
        # (per questo, per misurare una corrente, si può fare uso di generatori di tensione 
        # da 0V)
    
        nv_1 = len(circ.nodes_dict) - 1 # numero di soluzioni di tensione (al netto del ref)
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

        for elem in circ:
            if circuit.is_elem_voltage_defined(elem):
                index = index + 1
                varname = ("I("+elem.part_id.upper()+")").upper()
                self.variables += [varname]
                self.results.update({varname: x[index, 0]})
                self.errors.update({varname: error[index, 0]})
                self.units.update({varname: "A"})

        self.op_info = self.get_elements_op(circ, x)

    def __str__(self):
        str_repr = \
            "OP simulation results for %s (netlist %s).\nRun on %s, data filename %s.\n" % \
            (self.netlist_title, self.netlist_file, self.timestamp, self.filename)
        for v in self.variables:
            str_repr += "%s:\t%e\t%s\t(%e %s, %f %%)\n" % \
                (v, self.results[v], self.units[v], \
                self.errors[v], self.units[v], \
                self.errors[v]/self.results[v]*100.0)
        return str_repr

    def get_type(self):
        return "OP"

    def asmatrix(self):
        return self.x

    def get_table_array(self):
        table = [("Variable", "Value", "", "Error")]
        for v in self.variables:
            if self.results[v] != 0:
                relerror = self.errors[v]/self.results[v]*100.0
            else:
                relerror = 0.0
            line = (v, self.results[v], self.units[v],\
                "(% .2g %s, %.0f %%)" % (self.errors[v], 
                self.units[v], relerror))
            line = map(str, line)
            table.append(line)
        return table

    def get_elements_op(self, circ, x):
        tot_power = 0
        i_index = 0
        nv_1 = len(circ.nodes_dict) - 1
        op_info = []
        for elem in circ:
            ports_v_v = []
            if hasattr(elem, "get_op_info"):
                if elem.is_nonlinear:
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
                op_info += [elem.get_op_info(ports_v_v)]
            if isinstance(elem, devices.GISource):
                v = 0
                v = v + x[elem.n1-1] if elem.n1 != 0 else v
                v = v - x[elem.n2-1] if elem.n2 != 0 else v
                vs = 0
                vs = vs + x[elem.n1-1] if elem.n1 != 0 else vs
                vs = vs - x[elem.n2-1] if elem.n2 != 0 else vs
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
            elif circuit.is_elem_voltage_defined(elem):
                i_index = i_index + 1

        op_info.append("TOTAL POWER: %e W\n" % (tot_power,))
        return op_info

    def write_to_file(self, filename=None):
        if filename is None and self.filename is None:
            # maybe warn the user here?
            return
        if filename is None:
            filename = self.filename
        if filename != 'stdout':
            fp = open(filename+"info", "w")
        else:
            fp = sys.stdout
        fp.write(self.timestamp+"\n")
        fp.write("ahkab v. "+VERSION+u" (c) 2006-2013 Giuseppe Venturini\n\n")
        fp.write("Operating Point (OP) analysis\n\n")
        fp.write("Netlist: %s\nTitle: %s\n" % (self.netlist_file, self.netlist_title))
        fp.write("At %.2f K\n" % (self.temp,))
        fp.write("Options:\n\tvea = %e\n\tver = %f\n\tiea = %e\n\tier = %f\n\tgmin = %e\n" \
                 % (self.vea, self.ver, self.iea, self.ier, self.gmin))
        fp.write("\nConvergence reached in %d iterations.\n" % (self.iterations,))
        fp.write("\nResults:\n")
        vtable = self.get_table_array()
        fp.write(printing.table_setup(vtable))
        fp.write("\nELEMENTS OP INFORMATION:\n")
        for opi in self.op_info:    
            fp.write(opi)
            fp.write("-------------------\n")
        fp.flush()
        if filename != 'stdout':
            fp.close()
            solution._add_data(self, self.x)

    def print_short(self):
        str_repr = ""
        for v in self.variables:
            str_repr += "%s: %e %s,\t" % \
                (v, self.results[v], self.units[v])
        print str_repr[:-2]

    @staticmethod
    def gmin_check(op2, op1):
        """Checks the differences between two sets of OP results.
        (It is assumed that one set of results is calculated with Gmin, the other without.)
        op1, op2: the results vectors, interchangeable
        
        THIS METHOS IS UNBOUNDED.

        Returns:
        The list of the variables that did not pass the test. They are extracted from 
        the op_solution objects. 
        An empty list, if the check was passed.
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
                print "Unrecognized unit... Bug."
        return check_failed_vars

class ac_solution(solution):
    def __init__(self, circ, ostart, ostop, opoints, stype, op, outfile):
        """Holds a set of AC results.
            circ: the circuit instance of the simulated circuit
            ostart: the sweep starting angular frequency
            ostop: the sweep stopping frequency
            stype: the sweep type
            op: the linearization op used to compute the results
        """
        solution.__init__(self, circ, outfile)
        self.linearization_op = op
        self.stype = stype
        self.ostart, self.ostop, self.opoints = ostart, ostop, opoints
    
        nv_1 = len(circ.nodes_dict) - 1 # numero di soluzioni di tensione (al netto del ref)
        self.variables += ["w"]
        self.units.update({"w": "rad/s"})
    
        for index in range(nv_1):
            varname_abs = "|V%s|" % (str(circ.nodes_dict[index + 1]),)
            varname_arg = "arg(V%s)" % (str(circ.nodes_dict[index + 1]),)
            self.variables += [varname_abs]
            self.variables += [varname_arg]
            self.units.update({varname_abs: "V"})
            self.units.update({varname_arg: ""})
            if circ.is_int_node_internal_only(index+1):
                self.skip_nodes_list.append(index)

        for elem in circ: 
            if circuit.is_elem_voltage_defined(elem):
                varname_abs = "|I(%s)|" % (elem.part_id.upper(),)
                varname_arg = "arg(I(%s))" % (elem.part_id.upper(),)
                self.variables += [varname_abs]
                self.variables += [varname_arg]
                self.units.update({varname_abs: "A"})
                self.units.update({varname_arg: ""})

    def __str__(self):
        return "<AC simulation results for %s (netlist %s). %s sweep, from %g Hz to %g Hz, \
%d points. Run on %s, data filename %s.>" % \
        (self.netlist_title, self.netlist_file, self.stype, self.ostart, self.ostop, self.opoints, self.timestamp, self.filename)

    def add_line(self, omega, x):
        omega = numpy.mat(numpy.array([omega]))

        xsplit = numpy.mat(numpy.zeros((x.shape[0]*2, 1)))
        for i in range(x.shape[0]):
            xsplit[2*i, 0] = numpy.abs(x[i, 0])
            xsplit[2*i+1, 0] = numpy.angle(x[i, 0], deg=options.ac_phase_in_deg)
             
        data = numpy.concatenate((omega, xsplit), axis=0)
        solution._add_data(self, data)

    def get_type(self):
        return "AC"
        
    def get_x(self):
        return self.get(self.variables[0])

    def get_xlabel(self):
        return self.variables[0]

class dc_solution(solution):
    def __init__(self, circ, start, stop, sweepvar, stype, outfile):
        """Holds a set of DC results.
            circ: the circuit instance of the simulated circuit
            start: the sweep start value
            stop: the sweep stop value
            sweepvar: the swept variable
            stype: type of sweep
            outfile: the file to write the results to. 
                     (Use "stdout" to write to std output)
        """
        solution.__init__(self, circ, outfile)
        self.start, self.stop = start, stop
        self.stype = stype
    
        nv_1 = len(circ.nodes_dict) - 1 # numero di soluzioni di tensione (al netto del ref)
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
        return "<DC simulation results for %s (netlist %s). %s sweep of %s from %g %s to \
%g %s. Run on %s, data filename %s.>" % \
        (
         self.netlist_title, self.netlist_file, self.stype, self.variables[0].upper(), \
         self.start, self.units[self.variables[0]], self.stop, self.units[self.variables[0]], 
         self.timestamp, self.filename
        )

    def add_op(self, sweepvalue, op):
        """A DC sweep is made of a set of OP points. This method adds an OP solution and 
        its corresponding sweep value to the results set.
        """
        sweepvalue = numpy.mat(numpy.array([sweepvalue]))
        x = op.asmatrix()
        data = numpy.concatenate((sweepvalue, x), axis=0)
        solution._add_data(self, data)

    def get_type(self):
        return "DC"
        
    def get_x(self):
        return self.get(self.variables[0])

    def get_xlabel(self):
        return self.variables[0]

class tran_solution(solution):
    def __init__(self, circ, tstart, tstop, op, method, outfile):
        """Holds a set of TRANSIENT results.
            circ: the circuit instance of the simulated circuit
            tstart: the sweep starting angular frequency
            tstop: the sweep stopping frequency
            op: the op used to start the tran analysis
            outfile: the file to write the results to. 
                     (Use "stdout" to write to std output)
        """
        solution.__init__(self, circ, outfile)
        self.start_op = op
        self.tstart, self.tstop = tstart, tstop
        self.method = method

        self._lock = False
    
        nv_1 = len(circ.nodes_dict) - 1 # numero di soluzioni di tensione (al netto del ref)
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
        return "<TRAN simulation results for %s (netlist %s), from %g s to %g s. Diff. \
method %s. Run on %s, data filename %s.>" % \
        (
         self.netlist_title, self.netlist_file, self.tstart, self.tstop, self.method, 
         self.timestamp, self.filename
        )

    def add_line(self, time, x):
        """This method adds a solution and its corresponding time value to the results set.
        """
        if not self._lock:
            time = numpy.mat(numpy.array([time]))
            data = numpy.concatenate((time, x), axis=0)
            solution._add_data(self, data)
        else:
            printing.print_general_error(
                                "Attempting to add values to a complete result set. BUG"
                                )

    def lock(self):
        self._lock = True

    def get_type(self):
        return "TRAN"
        
    def get_x(self):
        return self.get(self.variables[0])

    def get_xlabel(self):
        return self.variables[0]


class pss_solution(solution):
    def __init__(self, circ, method, period, outfile, t_array=None, x_array=None):
        """Holds a set of TRANSIENT results.
            circ: the circuit instance of the simulated circuit
            tstart: the sweep starting angular frequency
            tstop: the sweep stopping frequency
            op: the op used to start the tran analysis
            outfile: the file to write the results to. 
                     (Use "stdout" to write to std output)
        """
        solution.__init__(self, circ, outfile)
        self.period = period
        self.method = method

        # We have mixed current and voltage results
        nv_1 = len(circ.nodes_dict) - 1 # numero di soluzioni di tensione (al netto del ref)
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

        if t_array is not None and x_array is not None:
            self.set_results(t_array, x_array)

    def __str__(self):
        return "<PSS simulation results for %s (netlist %s), period %g s. Method: %s. \
Run on %s, data filename %s.>" % \
        ( 
         self.netlist_title, self.netlist_file, self.period, self.method, self.timestamp, 
         self.filename
        )

    def set_results(self, t, x):
        """All the results are set at the same time for a PSS"""
        time = numpy.mat(numpy.array(t))
        data = numpy.concatenate((time, x), axis=0)
        solution._add_data(self, data)

    def asmatrix(self, verbose=3):
        allvalues = csvlib.load_csv(self.filename, load_headers=[], nsamples=None, skip=0L, verbose=verbose)
        return allvalues[0,:], allvalues[1:,:]

    def get_type(self):
        return "PSS"
        
    def get_x(self):
        return self.get(self.variables[0])

    def get_xlabel(self):
        return self.variables[0]

class symbolic_solution():
    def __init__(self, results_dict, substitutions, circ, outfile=None, tf=False):
        """Holds a set of Symbolic results.
            results_dict: the results dict returned by sympy.solve(),
            substitutions: the substitutions (dict) employed before solving,
            circ: the circuit instance of the simulated circuit.
            tf: is this set of results a set of transfer functions?
        """
        self.timestamp = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
        self.netlist_file = circ.filename
        self.netlist_title = circ.title
        self.substitutions = substitutions
        self.tf = tf

        # the keys are strings
        # self.symbols = map(str, results_dict.keys())
        self.results = case_insensitive_dict()
        for symbol, result in results_dict.iteritems():
            self.results.update({str(symbol):result})
        
        self._symbols = results_dict.keys() # keep them, they're useful
        for expr in results_dict.values():
            if tf:
                expr = expr['gain']
            for symb in expr.atoms():
                if symb.is_Symbol and not symb in self._symbols:
                    self._symbols.append(symb)
        self.filename = outfile if outfile != 'stdout' else None
        if self.filename is not None:
            self.save()

    def as_symbol(self, variable):
        symbs = filter(lambda x: x.name.lower() == variable.lower(), self._symbols)
        if len(symbs) == 0:
            raise ValueError, "No symbol %s in the results set."%(variable,)
        else:
            return symbs[0]
            
    def as_symbols(self, spacedstr):
        return map(self.as_symbol, spacedstr.split())
        
    def save(self):
        if not self.filename:
            raise Exception, "Writing the results to file requires setting the \
                              'filename' attribute"
        with open(self.filename, 'wb') as fp:
            pickle.dump(self, fp)
    
    @staticmethod
    def load(filename):
        with open(filename, 'rb') as fp:
            asolution = pickle.load(fp)
        return asolution


    def __repr__(self):
        return self.results.__repr__()

    def __str__(self):
        str_repr = "Symbolic %s results for %s (netlist %s).\nRun on %s.\n" % \
                   ('simulation'*(not self.tf) + 'transfer function'*self.tf, 
                    self.netlist_title, self.netlist_file, self.timestamp)
        keys = self.results.keys()
        keys.sort(lambda x, y: cmp(str(x), str(y))) 
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

    def get_type(self):
        return "Symbolic"

    # Access as a dictionary:
    def __len__(self):
        """Get the number of variables in the results set."""
        return len(self.results)

    def __getitem__(self, name):
        """Get a specific header, as from a dictionary."""
        return self.results[str(name).upper()]

    def get(self, name, default=None):
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
        return self.results.keys()

    def values(self):
        """Get all of the results set's variable's values."""
        return self.results.values()

    def items(self):
        return self.results.items()

    # iterator methods
    def __iter__(self):
        return self

    def next(self):
        if self.iter_index == len(self.results.keys())-1:
            self.iter_index = 0
            raise StopIteration
        else:
            self.iter_index += 1
        return self.results.keys()[self.iter_index], \
               self.results[self._symbols[self.iter_index]]


class case_insensitive_dict():
    """A dictionary that returns the same values for __str__.lower(key) and __str__.upper(key)
    """
    def __init__(self):
        self._dict = {}
        # Access as a dictionary BY VARIABLE NAME:
    def __len__(self):
        """Get the number of elements in the set."""
        return len(self._dict)

    def __repr__(self):
        rpr = "{"
        keys = self._dict.keys()
        for i in range(len(keys)):
            rpr += "%s: %s" % (keys[i], self._dict[keys[i]])
            if i < len(keys) - 1:
                rpr += ", "
            else:
                rpr += "}"
        return rpr

    def __getitem__(self, name):
        """Get a specific variable, as from a dictionary."""
        keys = self._dict.keys()
        i = map(str.upper, keys).index(name.upper())
        return self._dict[keys[i]]

    def get(self, name, default=None):
        """Given the key 'name', return its corresp. value. If not found, return 'default'
        """
        try:
            keys = self._dict.keys()
            i = map(str.upper, keys).index(name.upper())
        except KeyError:
            return default
        return self._dict[keys[i]]

    def has_key(self, name):
        """Determine whether the result set contains a variable."""
        return name.upper() in map(str.upper, self._dict.keys())

    def __contains__(self, name):
        """Determine whether the result set contains a variable."""
        return name.upper() in map(str.upper, self._dict.keys())

    def keys(self):
        """Get all keys"""
        return self._dict.keys()

    def values(self):
        """Get all values"""
        return self._dict.values()

    def items(self):
        """Get all keys and values pairs"""
        return self._dict.items()
        
    def update(self, adict):
        """Update the dictionary contents with the dictionary 'adict'"""
        return self._dict.update(adict)

    # iterator methods
    def __iter__(self):
        """Iterator"""
        return self._dict.__iter__()

