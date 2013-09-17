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

""" ahkab is an easy electronic circuit simulator.
"""

import sys
import tempfile
from optparse import OptionParser

import numpy
import sympy
import matplotlib

import ahkab

# analyses
import dc_analysis
import transient
import ac
import pss
import symbolic

import netlist_parser

import options
import constants 
import utilities

import plotting
import printing 

__version__ = "0.06a"

global _queue, _x0s, _print, _of

_queue = []
_print = False
_x0s = {None:None}
_of = []

def new_op(guess=True, x0=None, outfile=None, verbose=None):
	"""Assembles an OP analysis and returns the analysis object.
	The analysis itself can be run with:
	ahkab.run(...)
	or queued and then run subsequently.
	
	Parameters:
	guess (boolean): if set to True, the analysis will start from an initial guess,
	                 hopefully speeding up the convergence of stiff circuits.
	x0 (numpy matrix): if the 'guess' option above is not used, one can provide 
	                   a starting point directly, setting x0 to an opportunely sized
	                   numpy array. FIXME help method here
	                   If both x0 and guess are set, x0 takes the precedence.
	outfile (string): the filename of the output file where the results will be written.
	                  '.op' is automatically added at the end to prevent different 
	                  analyses from overwriting each-other's results.
	verbose (int): the verbosity level, from 0 (silent) to 6 (debug).
	
	Returns: the analysis object (a dict)
	"""
	if outfile is None or outfile == 'stdout': 
		if options.cli:
			outfile = 'stdout'
		else:
			ofi, outfile = tempfile.mkstemp(suffix='.op')
			_of.append(ofi) # keep the file open until quitting
	else:
		outfile += '.op'
	if verbose is None:
		verbose = 0
	return {'type':'op', 'guess':guess, 'x0':x0, 'outfile':outfile, 'verbose':verbose}
	
def new_dc(start, stop, points, source, sweep_type='LINEAR', guess=True, x0=None, outfile=None, \
           verbose=None):
	"""Assembles a DC sweep analysis and returns the analysis object.

	The analysis itself can be run with:
	ahkab.run(...)
	or queued and then run subsequently.
	
	Parameters:
	start (float): the start value for the sweep
	stop (float): the stop value for the sweep, included in the sweep
	points (int): the number of sweep points
	source (string): the independent current or voltage source to be swept. 
	sweep_type (string): can be set to either options.dc_lin_step (LINEAR) or 
	                     options.dc_log_step (LOG). Defaults to linear.
	guess (boolean): if set to True, the analysis will start from an initial guess,
	                 hopefully speeding up the convergence of stiff circuits.
	x0 (numpy matrix): if the 'guess' option above is not used, one can provide 
	                   a starting point directly, setting x0 to an opportunely sized
	                   numpy array. FIXME help method here
	                   If both x0 and guess are set, x0 takes the precedence.
	outfile (string): the filename of the output file where the results will be written.
	                  '.dc' is automatically added at the end to prevent different 
	                  analyses from overwriting each-other's results.
	verbose (int): the verbosity level, from 0 (silent) to 6 (debug).
	
	Returns: the analysis object (a dict)
	"""
	if outfile is None or outfile == 'stdout': 
		if options.cli:
			outfile = 'stdout'
		else:
			ofi, outfile = tempfile.mkstemp(suffix='.dc')
			_of.append(ofi) # keep the file open until quitting
	else:
		outfile += '.dc'
	if verbose is None:
		verbose = 0
	return {'type':'dc', 'start':float(start), 'stop':float(stop), 'points':float(points), 
	        'source':source, 'x0':x0, 'outfile':outfile, 'guess':guess, 'sweep_type':sweep_type, 
	        'verbose':verbose}
	
def new_tran(tstart, tstop, tstep, x0='op', method=transient.TRAP, use_step_control=True, 
             outfile=None, verbose=None):
	"""Assembles a TRAN analysis and returns the analysis object.

	The analysis itself can be run with:
	ahkab.run(...)
	or queued with ahakab.queue(...) and then run subsequently.
	
	Parameters:
	tstart (float): the start time for the transient analysis
	tstop (float): its stop time
	tstep (float): the time step. If the step control is active, this is the 
	               minimum time step value.
	x0 (numpy matrix): the optional starting point x0 = x(t=0). 
	method (string): the differentiation method to be used. Can be set to 
	                 'IMPLICIT_EULER', 'TRAP', 'GEAR4', 'GEAR5' or 'GEAR6'.
	                 Defaults to 'TRAP'.
	use_step_control (boolean): if False, use a fixed time step equal to tstep.
	outfile (string): the filename of the output file where the results will be written.
	                  '.tran' is automatically added at the end to prevent different 
	                  analyses from overwriting each-other's results.
	verbose (int): the verbosity level, from 0 (silent) to 6 (debug).
	
	Returns: the analysis object (a dict)
	"""
	if outfile is None or outfile == 'stdout': 
		if options.cli:
			outfile = 'stdout'
		else:
			ofi, outfile = tempfile.mkstemp(suffix='.tran')
			_of.append(ofi) # keep the file open until quitting
	else:
		outfile += '.tran'
	if verbose is None:
		verbose = 0
	return {"type":"tran", "tstart":tstart, "tstop":tstop, "tstep":tstep, 
	       "method":method, "use_step_control":use_step_control, 'x0':x0, 
	       'outfile':outfile, 'verbose':verbose}
	
def new_ac(start, stop, points, x0='op', sweep_type='LOG', outfile=None, verbose=None):
	"""Assembles an AC analysis and returns the analysis object.

	The analysis itself can be run with:
	ahkab.run(...)
	or queued with ahakab.queue(...) and then run subsequently.
	
	Parameters:
	start (float): the start angular frequency for the AC analysis
	stop (float): stop angular frequency
	points (float): the number of points to be use the discretize the 
					[start, stop] interval.
	sweep_type (string): Either 'LOG' or 'LINEAR', defaults to 'LOG'.
	outfile (string): the filename of the output file where the results will be written.
	                  '.ac' is automatically added at the end to prevent different 
	                  analyses from overwriting each-other's results.
	                  If unset or set to None, defaults to stdout.
	verbose (int): the verbosity level, from 0 (silent) to 6 (debug).
	
	Returns: the analysis object (a dict)
	"""
	if outfile is None or outfile == 'stdout': 
		if options.cli:
			outfile = 'stdout'
		else:
			ofi, outfile = tempfile.mkstemp(suffix='.ac')
			_of.append(ofi) # keep the file open until quitting
	else:
		outfile += '.ac'
	if verbose is None:
		verbose = 0
	return {'type':'ac', 'start':start, 'stop':stop, 'points':points, 'sweep_type':sweep_type, 
	        'x0':x0, 'outfile':outfile, 'verbose':verbose}
			
def new_pss(period, x0=None, points=None, method='brute-force', autonomous=False, outfile=None, 
            verbose=None):
	"""Assembles a Periodic Steady State (PSS) analysis and returns the analysis object.

	The analysis itself can be run with:
	ahkab.run(...)
	or queued with ahakab.queue(...) and then run subsequently.

	Parameters:
	period (float): the time period of the solution
	x0 (numpy matrix): the starting point solution, used at t=0.
	points (int): the number of points to use to discretize the PSS solution
	              If not set, if method is 'shooting', defaults to 
	              options.shooting_max_nr_iter
	method (string): can be either ahkab.BFPSS or ahkab.SHOOTING.
	autonomous (boolean): whether the circuit is autonomous or not.
	                      Non-autonomous circuits are currently unsupported!
	mna, Tf, D (numpy matrices): the numpy matrices to be used to solve the circuit.
	                             they are optional, but if they have already been 
	                             computed reusing them saves time.
	outfile (string): the filename of the output file where the results will be written.
	                  '.tran' is automatically added at the end to prevent different 
	                  analyses from overwriting each-other's results.
                      If unset defaults to stdout.
	verbose (int): the verbosity level, from 0 (silent) to 6 (debug).

	Returns: the analysis object (a dict)
	"""
	if outfile is None or outfile == 'stdout': 
		if options.cli:
			outfile = 'stdout'
		else:
			ofi, outfile = tempfile.mkstemp(suffix='.'+method.lower())
			_of.append(ofi) # keep the file open until quitting
	else:
		outfile += '.' + method.lower()
	if verbose is None:
		verbose = 0
	return {'type':"pss", "method":method, 'period':period, 'points':points,
	        'autonomous':autonomous, 'x0':x0, 'outfile':outfile, 'verbose':verbose}

def new_symbolic(source=None, ac_enable=True, r0s=False, subs=None, outfile=None, verbose=None):
	"""Assembles a Symbolic analysis and returns the analysis object.

	The analysis itself can be run with:
	ahkab.run(...)
	or queued with ahakab.queue(...) and then run subsequently.

	Parameters:
	source (string): if source is set, the transfer function between 'source'
	                 and each expression will be evaluated, including poles
	                 and zeros extraction. 'source' is to be set to the name
	                 of am independent current or voltage source present in 
	                 the circuit, eg. 'V1' or 'Iin'.
	ac_enable (boolean): if set True (default) the frequency-dependent elements will
	              be considered, otherwise the algorithm will focus on 
	              DC solutions (usually easier).
	r0s (boolean): if set to True, the finite output conductances of 
	               transistors go (where go = 1/r0) will be taken into 
	               account, otherwise it will be considered infinite (default).
	               In complex circuits, considering all gos may significanly 
	               complicate the solution - consider explicitly introducing
	               them where needed.
	subs (dict): is a dictionary of substitutions to be performed before 
	             attempting to solve the circuit. For example if two 
	             resistances R1 and R2 are to be equal, set subs={'R2':'R1'}
	             and R1 will be replaced by an instance of R2. This may 
	             simplify the solution (or allow finding one in reasonable 
	             time).
	outfile (string): the filename of the output file where the results will
	                  be written. '.symbolic' is automatically added at the 
	                  end to prevent different analyses from overwriting 
	                  each-other's results. If unset defaults to stdout.
	verbose (int): the verbosity level, from 0 (silent) to 6 (debug).

	Returns: the analysis object (a dict)
	"""
	if outfile is None or outfile == 'stdout': 
		if options.cli:
			outfile = 'stdout'
		else:
			ofi, outfile = tempfile.mkstemp(suffix='.symbolic')
			_of.append(ofi) # keep the file open until quitting
	else:
		outfile += '.symbolic'
	if verbose is None:
		verbose = 0
	return {'type':"symbolic", 'source':source, 'ac_enable':ac_enable, 'r0s':r0s, 'subs':subs, 
		'outfile':outfile, 'verbose':verbose}
		
def queue(*analysis):
	global _queue
	for an in analysis: # let's hope the user knows what he's doing!
		_queue += [an]
	
def run(circ, an_list=None):
	""" Processes an analysis vector:
	circ: the circuit instance
	an_queue: the list of analyses to be performed, if unset defaults to those queued
		      with queue_analysis()
	
	Returns: the results (in dict form)
	"""
	results = {}
	
	if not an_list:
		an_list = _queue
	
	while len(an_list):
		an_item = an_list.pop(0)
		an_type = an_item.pop('type')
		r = analysis[an_type](circ, **an_item)
		results.update({an_type:r})
		if an_type == 'op':
			_x0s.update({'op':r})
			_x0s.update({'op+ic':icmodified_x0(circ, r)})
			_handle_netlist_ics(circ, an_list, ic_list=[])
	return results

def new_x0(circ, icdict):
	return dc_analysis.build_x0_from_user_supplied_ic(circ, icdict)

def icmodified_x0(circ, x0):
	return dc_analysis.modify_x0_for_ic(circ, x0)
	
def get_op_x0(circ):
	return run(circ, [new_op()])

def set_temperature(T):
	T = float(T)
	if T > 300: 
		printing.print_warning(u"The temperature will be set to %f \xB0 C.")
	constants.T = utilities.Celsius2Kelvin(T)

def process_postproc(postproc_list, title, results, outfilename):
	"""Runs the post-processing operations, such as plotting.
	postproc_list: list of post processing operations as returned by main()
	title: the deck title
	results: the results to be plotted (including the ones that are not needed)
	outfilename: if the plots are saved to disk, this is the filename without extension

	Returns: None
	"""
	index = 0
	if outfilename == 'stdout':
		printing.print_warning("Plotting and printing the results to stdout are incompatible options. Plotting skipped.")
		return
	for postproc in postproc_list:
		plotting.plot_results(title, postproc["x"], postproc["l2l1"], results[postproc["analysis"]], "%s-%d.%s" % (outfilename, index, options.plotting_outtype))
		index = index +1
	if len(postproc_list) and options.plotting_show_plots:
		plotting.show_plots()
	return None
	
analysis = {'op':dc_analysis.op_analysis, 'dc': dc_analysis.dc_analysis, 
            'tran': transient.transient_analysis, 'ac': ac.ac_analysis,
            'pss': pss.pss_analysis, 'symbolic': symbolic.symbolic_analysis, 
            'temp':set_temperature}

def main(filename, outfile="stdout", verbose=3):
	"""This method allows calling ahkab from a Python script.
	"""
	printing.print_info_line(("This is ahkab %s running with:" %(__version__), 6), verbose)
	printing.print_info_line(("  Python %s" % (sys.version.split('\n')[0],), 6), verbose)
	printing.print_info_line(("  Numpy %s"  % (numpy.__version__), 6), verbose)
	printing.print_info_line(("  Sympy %s"  % (sympy.__version__), 6), verbose)
	printing.print_info_line(("  Matplotlib %s"  % (matplotlib.__version__), 6), verbose)

	read_netlist_from_stdin = (filename is None or filename == "-")
	(circ, directives, postproc_direct) = netlist_parser.parse_circuit(filename, read_netlist_from_stdin)
	
	check, reason = dc_analysis.check_circuit(circ)
	if not check:
		printing.print_general_error(reason)
		printing.print_circuit(circ)
		sys.exit(3)
	
	if verbose > 3 or _print:
		print "Parsed circuit:"
		printing.print_circuit(circ)

	ic_list = netlist_parser.parse_ics(directives)
	_handle_netlist_ics(circ, an_list=[], ic_list=ic_list)
	results = {}
	for an in netlist_parser.parse_analysis(circ, directives):
		if 'outfile' not in an.keys() or not an['outfile']:
			an.update({'outfile':outfile+("."+an['type'])*(outfile != 'stdout')})
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
		ic_label = ic.keys()[0]
		icdict = ic[ic_label]
		_x0s.update({ic_label:new_x0(circ, icdict)})
	for an in an_list:
		if isinstance(an['x0'], str):
			if an['x0'] in _x0s.keys():
				an['x0'] = _x0s[an['x0']]
			elif an_list.index(an) == 0:
				printing.print_general_error("Unknown x0 %s" % an["x0"])
				sys.exit(54)


if __name__ == "__main__":
	usage = "usage: \t%prog [options] <filename>\n\nThe filename is the netlist to be open. Use - (a dash) to read from stdin."
	parser = OptionParser(usage, version="%prog "+__version__+" (c) 2006-2013 Giuseppe Venturini")
	
	#general options
	parser.add_option("-v", "--verbose", action="store", type="string", dest="verbose", default="3", help="Verbose level: from 0 (almost silent) to 5 (debug)")
	parser.add_option("-p", "--print", action="store_true", dest="print_circuit", default=False, help="Print the parsed circuit")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", default=options.default_output_file, help="Data output file. Defaults to stdout.")
	parser.add_option("", "--dc-guess", action="store", type="string", dest="dc_guess", default="guess", help="Guess to be used to start a op or dc analysis: none or guess. Defaults to guess.")
	parser.add_option("-t", "--tran-method", action="store", type="string", dest="method", default=transient.TRAP.lower(), help="Method to be used in transient analysis: " +transient.IMPLICIT_EULER.lower()+", "+transient.TRAP.lower()+", "+transient.GEAR2.lower()+", "+transient.GEAR3.lower()+", "+transient.GEAR4.lower()+", "+transient.GEAR5.lower()+" or "+transient.GEAR6.lower()+". Defaults to TRAP.")
	parser.add_option("", "--t-fixed-step", action="store_true", dest="no_step_control", default=False, help="Disables the step control in transient analysis. Useful if you want to perform a FFT on the results.")
	parser.add_option("", "--v-absolute-tolerance", action="store", type="string", dest="vea", default=None, help="Voltage absolute tolerance. Default: "+str(options.vea)+" V")
	parser.add_option("", "--v-relative-tolerance", action="store", type="string", dest="ver", default=None, help="Voltage relative tolerance. Default: "+str(options.ver))
	parser.add_option("", "--i-absolute-tolerance", action="store", type="string", dest="iea", default=None, help="Current absolute tolerance. Default: "+str(options.iea)+" A")
	parser.add_option("", "--i-relative-tolerance", action="store", type="string", dest="ier", default=None, help="Current relative tolerance. Default: "+str(options.ier))
	parser.add_option("", "--h-min", action="store", type="string", dest="hmin", default=None, help="Minimum time step. Default: "+str(options.hmin))
	parser.add_option("", "--dc-max-nr", action="store", type="string", dest="dc_max_nr_iter", default=None, help="Maximum nr of NR iterations for dc analysis. Default: "+str(options.dc_max_nr_iter))
	parser.add_option("", "--t-max-nr", action="store", type="string", dest="transient_max_nr_iter", default=None, help="Maximum nr of NR iterations for each time step during transient analysis. Default: "+str(options.transient_max_nr_iter))
	parser.add_option("", "--t-max-time", action="store", type="string", dest="transient_max_time_iter", default=None, help="Maximum nr of time iterations during transient analysis. Setting it to 0 (zero) disables the limit. Default: "+str(options.transient_max_time_iter))
	parser.add_option("", "--s-max-nr", action="store", type="string", dest="shooting_max_nr_iter", default=None, help="Maximum nr of NR iterations during shooting analysis. Setting it to 0 (zero) disables the limit. Default: "+str(options.shooting_max_nr_iter))
	parser.add_option("", "--gmin", action="store", type="string", dest="gmin", default=None, help="The minimum conductance to ground. Inserted when requested. Default: "+str(options.gmin))
	parser.add_option("", "--cmin", action="store", type="string", dest="cmin", default=None, help="The minimum capacitance to ground. Default: "+str(options.cmin))
	parser.add_option("", "--eps", action="store_true", dest="eps", default=False, help="Calculate the machine precision. The machine precision defaults to "+str(utilities.EPS))
	
	(cli_options, remaning_args) = parser.parse_args()
	
	verbose = int(cli_options.verbose)
	if cli_options.method is not None:
		options.default_tran_method = cli_options.method.upper()
	if cli_options.vea is not None:
		options.vea = float(cli_options.vea)
	if cli_options.ver is not None:
		options.ver = float(cli_options.ver)
	if cli_options.iea is not None:
		options.iea = float(cli_options.iea)
	if cli_options.ier is not None:
		options.ier = float(cli_options.ier)
	if cli_options.hmin is not None:
		options.hmin = float(cli_options.hmin)
	if cli_options.dc_max_nr_iter is not None:
		options.dc_max_nr_iter = int(float(cli_options.dc_max_nr_iter))
	if cli_options.transient_max_nr_iter is not None:
		options.transient_max_nr_iter = int(float(cli_options.transient_max_nr_iter))
	if cli_options.transient_max_time_iter is not None:
		options.transient_max_time_iter = int(float(cli_options.transient_max_time_iter))
	if cli_options.shooting_max_nr_iter is not None:
		options.shooting_max_nr_iter = int(float(cli_options.shooting_max_nr_iter))
	if cli_options.gmin is not None:
		options.gmin = float(cli_options.gmin)
	if cli_options.cmin is not None:
		options.cmin = float(cli_options.cmin)
	if cli_options.eps:
		utilities.EPS = utilities.calc_eps()
		print "Detected machine precision: " + str(utilities.EPS)

	if not len(remaning_args) == 1:
		print "Usage: ./ahkab.py [options] <filename>\n./ahkab.py -h for help"
		sys.exit(1)
	if remaning_args[0] == '-':
		read_netlist_from_stdin = True
	else:
		read_netlist_from_stdin = False
	if not read_netlist_from_stdin and not utilities.check_file(remaning_args[0]):
		sys.exit(23)

	options.transient_no_step_control = cli_options.no_step_control
	options.dc_use_guess = cli_options.dc_guess
	options.cli = True
	_print = cli_options.print_circuit

	# Program execution
	main(filename=remaning_args[0], outfile=cli_options.outfile, verbose=verbose)

	sys.exit(0)
