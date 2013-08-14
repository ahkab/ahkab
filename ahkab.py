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

""" ahkab is a easy electronic circuit simulator.
"""

import sys
from optparse import OptionParser

# analyses
import dc_analysis
import transient
import ac
import shooting
import bfpss
import symbolic

import netlist_parser

import options
import constants 
import utilities

import plotting
import printing 

VERSION = "0.06"


def process_analysis(an_list, circ, outfile, verbose, cli_tran_method=None, guess=True, disable_step_control=False):
	""" Processes an analysis vector:
	an_list: the list of analysis to be performed, as returned by netlist_parser
	circ: the circuit instance, returned by netlist_parser
	outfile: a filename. Results will be written to it. If set to stdout, prints to stdout
	verbose: verbosity level
	cli_tran_method: force the specified method in each tran analysis (see transient.py)
	guess: use the builtin method get_dc_guess to guess x0
	
	Returns: None
	"""
	x0_op = None
	x0_ic_dict = {}
	results = {}

	for directive in [ x for x in an_list if x["type"] == "ic" ]:
		x0_ic_dict.update({
			directive["name"]:\
			dc_analysis.build_x0_from_user_supplied_ic(circ, voltages_dict=directive["vdict"], currents_dict=directive["cdict"])
			})
	
	for an in an_list:
		if outfile != 'stdout':
			data_filename = outfile + "." + an["type"]
		else:
			data_filename = outfile

		if an["type"] == "ic":
			continue

		if an["type"] == "op":
			if an["guess_label"] is None:
				x0_op = dc_analysis.op_analysis(circ, guess=guess, data_filename=data_filename, verbose=verbose)
			else:
				if not an["guess_label"] in x0_ic_dict:
					printing.print_warning("op: guess is set but no matching .ic directive was found.")
					printing.print_warning("op: using built-in guess method: "+str(guess))
					x0_op = dc_analysis.op_analysis(circ, guess=guess, verbose=verbose)
				else:
					x0_op = dc_analysis.op_analysis(circ, guess=False, x0=x0_ic_dict[an["guess_label"]], verbose=verbose)
			sol = x0_op
		
		elif an["type"] == "dc":
			if an["source_name"][0].lower() == "v":
				elem_type = "vsource"
			elif an["source_name"][0].lower() == "i":
				elem_type = "isource"
			else:
				printing.print_general_error("Type of sweep source is unknown: " + an[1][0])
				sys.exit(1)
			sol = dc_analysis.dc_analysis(
					circ, start=an["start"], stop=an["stop"], step=an["step"], \
					type_descr=(elem_type, an["source_name"][1:]), 
					xguess=x0_op, data_filename=data_filename, guess=guess, 
					stype=an['stype'], verbose=verbose)
			
		
		#{"type":"tran", "tstart":tstart, "tstop":tstop, "tstep":tstep, "uic":uic, "method":method, "ic_label":ic_label}
		elif an["type"] == "tran":
			if cli_tran_method is not None:
				tran_method = cli_tran_method.upper()
			elif an["method"] is not None:
				tran_method = an["method"].upper()
			else:
				tran_method = options.default_tran_method
			
			# setup the initial condition (t=0) according to uic
			# uic = 0 -> all node voltages and currents are zero
			# uic = 1 -> node voltages and currents are those computed in the last OP analysis
			# uic = 2 -> node voltages and currents are those computed in the last OP analysis
			#            combined with the ic=XX directive found in capacitors and inductors
			# uic = 3 -> use a .ic directive defined by the user
			uic = an["uic"]
			if uic == 0:
				x0 = None
			elif uic == 1:
				if x0_op is None:
					printing.print_general_error("uic is set to 1, but no op has been calculated yet.")
					sys.exit(51)
				x0 = x0_op
			elif uic == 2:
				if x0_op is None:
					printing.print_general_error("uic is set to 2, but no op has been calculated yet.")
					sys.exit(51)
				x0 = dc_analysis.modify_x0_for_ic(circ, x0_op)
			elif uic == 3:
				if an["ic_label"] is None:
					printing.print_general_error("uic is set to 3, but param ic=<ic_label> was not defined.")
					sys.exit(53)
				elif not an["ic_label"] in x0_ic_dict:
					printing.print_general_error("uic is set to 3, but no .ic directive named %s was found." \
						%(str(an["ic_label"]),))
					sys.exit(54)
				x0 = x0_ic_dict[an["ic_label"]]
			
			sol = transient.transient_analysis(circ, \
				tstart=an["tstart"], tstep=an["tstep"], tstop=an["tstop"], \
				x0=x0, mna=None, N=None, verbose=verbose, data_filename=data_filename, \
				use_step_control=(not disable_step_control), method=tran_method)
		
		elif an["type"] == "shooting":
			if an["method"]=="brute-force":
				sol = bfpss.bfpss(circ, period=an["period"], step=an["step"], mna=None, Tf=None, \
					D=None, points=an["points"], autonomous=an["autonomous"], x0=x0_op, \
					data_filename=data_filename, verbose=verbose)
			elif an["method"]=="shooting":	
				sol = shooting.shooting(circ, period=an["period"], step=an["step"], mna=None, \
					Tf=None, D=None, points=an["points"], autonomous=an["autonomous"], \
					data_filename=data_filename, verbose=verbose)
		elif an["type"] == "symbolic":
			if not 'subs' in an.keys():
				an.update({'subs':None})
			sol = symbolic.solve(circ, an['source'], opts={'ac':an['ac']}, subs=an['subs'], verbose=verbose)
		elif an["type"] == "ac":
			sol = ac.ac_analysis(circ=circ, start=an['start'], nsteps=an['nsteps'], \
				stop=an['stop'], step_type='LOG', xop=x0_op, mna=None,\
			        data_filename=data_filename, verbose=verbose)
		elif an["type"] == "temp":
			constants.T = utilities.Celsius2Kelvin(an['temp'])
		results.update({an["type"]:sol})
	return results

def process_postproc(postproc_list, title, results, outfilename, remote=False):
	"""Runs the post-processing operations, such as plotting.
	postproc_list: list of post processing operations as returned by main()
	title: the deck title
	results: the results to be plotted (including the ones that are not needed)
	outfilename: if the plots are saved to disk, this is the filename without extension
	remote: boolean, do not show plots if True (such as ssh without X11 forwarding)

	Returns: None
	"""
	index = 0
	if outfilename == 'stdout':
		printing.print_warning("Plotting and printing the results to stdout are incompatible options. Plotting skipped.")
		return
	for postproc in postproc_list:
		#print postproc["analysis"], results.keys(), results.has_key(postproc["analysis"]), results[postproc["analysis"]] is None #DEBUG
		plotting.plot_results(title, postproc["x"], postproc["l2l1"], results[postproc["analysis"]], "%s-%d.%s" % (outfilename, index, options.plotting_outtype))
		index = index +1
	if len(postproc_list) and not remote:
		plotting.show_plots()
	return None

def main(filename, outfile="stdout", tran_method=transient.TRAP.lower(), no_step_control=False, dc_guess='guess', print_circuit=False, remote=True, verbose=3):
	"""This method allows calling ahkab from a Python script.
	"""

	utilities._set_execution_lock()

	read_netlist_from_stdin = (filename is None or filename == "-")
	(circ, directives, postproc_direct) = netlist_parser.parse_circuit(filename, read_netlist_from_stdin)
	
	check, reason = dc_analysis.check_circuit(circ)
	if not check:
		printing.print_general_error(reason)
		printing.print_circuit(circ)
		sys.exit(3)
	
	if verbose > 3 or print_circuit:
		print "Parsed circuit:"
		printing.print_circuit(circ)
	elif verbose > 1:
		print circ.title.upper()
	
	an_list = netlist_parser.parse_analysis(circ, directives)
	postproc_list = netlist_parser.parse_postproc(circ, an_list, postproc_direct)
	if len(an_list) > 0: 
		printing.print_info_line(("Requested an.:", 4), verbose)
		if verbose >= 4:
			map(printing.print_analysis, an_list)
	else:
		if verbose:
			printing.print_warning("No analysis requested.")
	if len(an_list) > 0:
		results = process_analysis(an_list, circ, outfile, verbose, guess=dc_guess.lower()=="guess", \
		cli_tran_method=tran_method, disable_step_control=no_step_control)
	else:
		printing.print_warning("Nothing to do. Quitting.")

	if len(an_list) > 0 and len(postproc_list) > 0 and len(results):
		process_postproc(postproc_list, circ.title, results, outfile, remote)

	utilities._unset_execution_lock()

	return results

if __name__ == "__main__":
	parser = OptionParser(usage="usage: \t%prog [options] <filename>\n\nThe filename is the netlist to be open. Use - (a dash) to read from stdin.",  version="%prog "+VERSION+u" (c) 2006-2013 Giuseppe Venturini")
	
	#general options
	parser.add_option("-v", "--verbose", action="store", type="string", dest="verbose", default="3", help="Verbose level: from 0 (almost silent) to 5 (debug)")
	parser.add_option("-p", "--print", action="store_true", dest="print_circuit", default=False, help="Print the parsed circuit")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", default="stdout", help="Data output file. Defaults to stdout.")
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
		method = cli_options.method.upper()
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
	
	# Program execution
	main(filename=remaning_args[0], outfile=cli_options.outfile, tran_method=cli_options.method, \
		no_step_control=cli_options.no_step_control, dc_guess=cli_options.dc_guess, \
		print_circuit=cli_options.print_circuit, remote=False, verbose=verbose)

	sys.exit(0)
