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
import netlist_parser, dc_analysis, transient, utilities, shooting, options, printing


VERSION = "0.02"


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
	last_x_tran = None

	for directive in [ x for x in an_list if x[0] == "ic" ]:
		x0_ic_dict.update({directive[1]:dc_analysis.build_x0_from_user_supplied_ic(circ,  voltages_dict=directive[2], currents_dict=directive[3])})
	
	for an in an_list:
		if outfile != 'stdout':
			data_filename = outfile + "." + an[0]
		else:
			data_filename = outfile

		if an[0] == "op":
			guess_label = an[1]
			if guess_label is None:
				x0_op = dc_analysis.op_analysis(circ, guess=guess, verbose=verbose)
			else:
				if not guess_label in x0_ic_dict:
					printing.print_warning("op: guess is set but no matching .ic directive was found.")
					printing.print_warning("op: using built-in guess method: "+str(guess))
					x0_op = dc_analysis.op_analysis(circ, guess=guess, verbose=verbose)
				else:
					x0_op = dc_analysis.op_analysis(circ, guess=False, x0=x0_ic_dict[guess_label], verbose=verbose)
		
		elif an[0] == "dc":
			if an[1][0].lower() == "v":
				elem_type = "vsource"
			elif an[1][0].lower() == "i":
				elem_type = "isource"
			else:
				printing.print_general_error("Type of sweep source is unknown: " + an[1][0])
				sys.exit(1)
			dc_analysis.dc_analysis(circ, start=an[2], stop=an[3], step=an[4], elem_type=elem_type, elem_descr=an[1][1:], data_filename=data_filename, guess=guess, verbose=verbose)
		
		elif an[0] == "tran":
			if cli_tran_method is not None:
				tran_method = cli_tran_method.upper()
			elif an[5] is not None:
				tran_method = an[5].upper()
			else:
				tran_method = options.default_tran_method
			
			# setup the initial condition according to uic
			# uic = 0 -> parti con tutto zero
			# uic = 1 -> parti con l'op, se disponibile.
			# uic = 2 -> parti con l'op, tieni conto delle ic di condensatori e induttori
			# uic = 3 -> carica un ic definito dall'utente
			uic = an[4]
			if uic == 0:
				x0 = None
			elif uic == 1:
				x0 = x0_op # if there's no x0_op, it's the same as uic=0
			elif uic == 2:
				if x0_op is None:
					printing.print_general_error("uic is set to 2, but no op has been calculated yet.")
					sys.exit(51)
				x0 = dc_analysis.modify_x0_for_ic(circ, x0_op)
			elif uic == 3:
				if an[6] is None:
					printing.print_general_error("uic is set to 3, but param ic=<ic_label> was not defined.")
					sys.exit(53)
				elif not an[6] in x0_ic_dict:
					printing.print_general_error("uic is set to 3, but no .ic directive named "+str(an[6])+" was found.")
					sys.exit(54)
				x0 = x0_ic_dict[an[6]]
			
			(time, last_x_tran) = transient.transient_analysis(circ, tstart=an[1], tstep=an[3], tstop=an[2], x0=x0, mna=None, N=None, verbose=verbose, data_filename=data_filename, use_step_control=(not disable_step_control), method=tran_method)
		
		elif an[0] == "shooting":
			shooting.shooting(circ, period=an[1], step=an[3], mna=None, Tf=None, D=None, points=an[2], autonomous=an[4], x0=last_x_tran, data_filename=data_filename, verbose=verbose)
				
	return None

if __name__ == "__main__":
	# Import Psyco if available
	try:
		import psyco
		psyco.full()
	except ImportError:
		pass
	
	parser = OptionParser(usage="usage: \t%prog [options] <filename>\n\nThe filename is the netlist to be open. Use - (a dash) to read from stdin.",  version="%prog "+VERSION+u" (c) 2006-2007 Giuseppe Venturini")
	
	#general options
	parser.add_option("-v", "--verbose", action="store", type="string", dest="verbose", default="3", help="Verbose level: from 0 (almost silent) to 5 (debug)")
	parser.add_option("-p", "--print", action="store_true", dest="print_circuit", default=False, help="Print the parsed circuit")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", default="stdout", help="Data output file. Defaults to stdout.")
	parser.add_option("", "--dc-guess", action="store", type="string", dest="dc_guess", default="guess", help="Guess to be used to start a op or dc analysis: none or guess. Defaults to guess.")
	parser.add_option("-t", "--tran-method", action="store", type="string", dest="method", default=None, help="Method to be used in transient analysis: " +transient.IMPLICIT_EULER.lower()+", "+transient.TRAP.lower()+", "+transient.GEAR2.lower()+", "+transient.GEAR3.lower()+", "+transient.GEAR4.lower()+", "+transient.GEAR5.lower()+" or "+transient.GEAR6.lower()+". Defaults to IE.")
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
	parser.add_option("", "--gmin", action="store", type="string", dest="gmin", default=None, help="The minimum conductance to ground to be inserted (when requested). Default: "+str(options.gmin))
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
	
	(circ, directives) = netlist_parser.parse_circuit(remaning_args[0], read_netlist_from_stdin)
	
	check, reason = dc_analysis.check_circuit(circ)
	if not check:
		printing.print_general_error(reason)
		printing.print_circuit(circ)
		sys.exit(3)
	
	if verbose > 3 or cli_options.print_circuit:
		print "Parsed circuit:"
		printing.print_circuit(circ)
	elif verbose > 1:
		print circ.title.upper()
	
	an_list = netlist_parser.parse_analysis(circ, directives)
	if verbose > 3:
		if len(an_list) > 0:
			print "Requested an.:"
			map(printing.print_analysis, an_list)
		else:
			printing.print_warning("No analysis requested.")
	if len(an_list) > 0:
		process_analysis(an_list, circ, cli_options.outfile, verbose, guess=cli_options.dc_guess.lower()=="guess", cli_tran_method=cli_options.method, \
		disable_step_control=cli_options.no_step_control)
	else:
		if verbose:
			printing.print_warning("Nothing to do. Quitting.")
	sys.exit(0)
