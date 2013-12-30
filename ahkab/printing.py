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
This is the printing module of the simulator. Using its functions, the output will
be somewhat uniform.
"""

import sys
import devices, options
import diode, mosq, ekv, switch
    
def print_circuit(circ):
    """Prints the whole circuit to stdout, in netlist format.
    
    Parameters:
    circ: the circuit instance to be printed.
    
    Returns: None
    """
    if circ.title:
        print circ.title
        
    for elem in circ:
        print_netlist_elem_line(elem, circ)
    
    return None
    
def print_netlist_elem_line(elem, circ):
    """Prints a elem to stdout from the provided circuit instance.
    
    Parameters:
    elem: the elem to be printed
    circ: the circuit instance to which the element belongs.
    
    Returns: None
    """
    if hasattr(elem, "n1") and hasattr(elem, "n2"):
        ext_n1 = circ.nodes_dict[elem.n1]
        ext_n2 = circ.nodes_dict[elem.n2]
    
    sys.stdout.write(elem.letter_id.upper() + elem.descr + " ")
    
    if isinstance(elem, devices.Resistor)   \
    or isinstance(elem, diode.diode)        \
    or isinstance(elem, devices.ISource)    \
    or isinstance(elem, devices.VSource)    \
    or isinstance(elem, devices.Capacitor)  \
    or isinstance(elem, devices.Inductor):
        sys.stdout.write(ext_n1 + " " + ext_n2 + " ")
    
    elif isinstance(elem, devices.EVSource) \
    or isinstance(elem, devices.GISource)   \
    or isinstance(elem, switch.switch_device):
        sys.stdout.write(ext_n1 + " " + ext_n2 + " " + circ.nodes_dict[elem.sn1]+ " " + \
        circ.nodes_dict[elem.sn2] + " ")

    elif isinstance(elem, devices.InductorCoupling):
        sys.stdout.write(" ")

    elif isinstance(elem, mosq.mosq_device): #  quadratic mos
        sys.stdout.write(ext_n1 + " " + circ.nodes_dict[elem.ng] + " " + ext_n2 + " ")
    
    elif isinstance(elem, ekv.ekv_device):
        sys.stdout.write(ext_n1 + " " + circ.nodes_dict[elem.ng] + " " + ext_n2 + " " + circ.nodes_dict[elem.nb] + " ")
    
    elif elem.letter_id == "y":
        sys.stdout.write(ext_n1 + " " + ext_n2 + " ")
    
    else:
        print_general_error("\nUnknown element, this is probably a bug: " + elem.__class__.__name__)
        sys.exit(1)
    
    #   Defaults to `Component.__str__`
    print str(elem)
    
    return None
    
def print_analysis(an):
    """Prints a analysis to stdout, with the netlist syntax
    
    Parameters:
    an: an analisys, a element of the list returned from netlist_parser.parse_analysis
    
    Returns: None
    """
    if an["type"] == "op":
        print ".op"
    elif an["type"] == "dc":
        print ".dc", an["source_name"], "start =", an["start"], "stop =", an["stop"], "step =", an["step"], "type =", an["stype"]
    elif an["type"] == "tran":
        sys.stdout.write(".tran tstep="+str(an["tstep"])+" tstop="+str(an["tstop"])+" tstart="+str(an["tstart"]))
        if an["method"] is not None:
            print " method=" + an["method"]
        else:
            print ""
    elif an["type"] == "shooting":
        sys.stdout.write(".shooting period="+ str(an["period"])+" method="+str(an["method"]))
        if an["points"] is not None:
            sys.stdout.write(" points=" + str(an["points"]))
        if an["step"] is not None:
            sys.stdout.write(" step=" + str(an["step"]))
        print " autonomous=", an["autonomous"]

def print_general_error(description, print_to_stdout=False):
    """Prints a error message to stderr.
    
    Parameters:
    description: the error's description
    print_to_stdout:
    
    Returns: None
    """
    the_error_message = "E: " + description
    if print_to_stdout:
        print the_error_message
    else:
        sys.stderr.write(the_error_message+"\n")
    return None

def print_warning(description, print_to_stdout=False):
    """Prints a warning message to stderr.
    
    Parameters:
    description: the warning's description
    print_to_stdout:
    
    Returns: None
    """
    the_warning_message = "W: " + description
    if print_to_stdout:
        print the_warning_message
    else:
        sys.stderr.write(the_warning_message+"\n")
    return None
    
def print_info_line((msg, relevance), verbose, print_nl=True):
    if verbose >= relevance:
        if print_nl:
            print msg
        else:
            print msg,
    # suppressed.

def print_parse_error(nline, line, print_to_stdout=False):
    """Prints a parsing error in the netlist to stderr.
    
    Parameters:
    nline: number of the line on which the error was found
    line: the line of the file
    print_to_stdout:
    
    Returns: None
    """
    print_general_error("Parse error on line " + str(nline) + ":", print_to_stdout)
    if print_to_stdout:
        print line
    else:
        sys.stderr.write(line+"\n")
    return None

def print_symbolic_results(x):
    keys = x.keys()
    keys.sort() 
    for key in keys:
        print str(key) + "\t = " + str(x[key])
    return None

def print_symbolic_transfer_functions(x):
    keys = x.keys()
    keys.sort() 
    for key in keys:
        print str(key) + " = " + str(x[key]['gain'])
        print '\tDC: ' + str(x[key]['gain0'])
        for index in range(len(x[key]['poles'])):
            print '\tP'+str(index)+":", str(x[key]['poles'][index])
        for index in range(len(x[key]['zeros'])):
            print '\tZ'+str(index)+":", str(x[key]['zeros'][index])
    return None
def print_symbolic_equations(eq_list):
    print "+--" 
    for eq in eq_list:
        print "| " + str(eq)
    print "+--"
    return

def print_result_check(badvars, verbose=2):
    """Prints out the results of the OP check performed by results.op_solution.gmin_check
    It assumes one set of results is calculated with Gmin, the other without.
    badvars: the list returned by results.op_solution.gmin_check
    
    Returns: None
    """
    if len(badvars):
        print "Warning: solution is heavvily dependent on gmin."
        print "Affected variables:"
        for bv in badvars:
            print bv
    else:
        if verbose: 
            print "Difference check is within margins." 
            print "(Voltage: er=" + str(options.ver) + ", ea=" + str(options.vea) + \
            ", Current: er=" + str(options.ier) + ", ea=" + str(options.iea) + ")"
    return None

def table_print(twodarray, separator='  '):
    print table_setup(twodarray, separator=separator)

def table_setup(twodarray, separator='  '):
    table_string = ""
    col_width = []
    if len(twodarray) == 0 or len(twodarray[0]) == 0:
        return
    for ci in range(len(twodarray[0])):
        current_width = 0
        for ri in range(len(twodarray)):
            elem_width = len(str(twodarray[ri][ci]))            
            if elem_width > current_width:
                current_width = elem_width
        col_width.append(current_width)
    for ri in range(len(twodarray)):
        current_str = "" 
        for ci in range(len(twodarray[ri])):
            elem = str(twodarray[ri][ci])
            elem_width = len(elem)
            if not ci +1 % 3 == 1:
                current_str = current_str + " "*(col_width[ci]-elem_width) + elem + separator
            else:
                current_str = current_str + elem + " "*(col_width[ci]-elem_width) +  separator
        table_string += current_str + "\n"
    return table_string

