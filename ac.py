# -*- coding: iso-8859-1 -*-
# ac.py
# AC analysis module
# Copyright 2010 Giuseppe Venturini
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

""" This module offers the methods required to perform an AC analysis.
Our problem can be written as:
	MNA*x + AC*x + J*x + Nac = 0
We need:
	1. the mna matrix MNA
	2. the AC matrix, holding the frequency dep parts
	3. The linearized non-linear elements end up in J
	4. Nac, the AC sources contribution  
In order for an AC analysis to be performed, an OP has to be computed first,
if there is any non-linear device in the circuit.
"""

import sys
import numpy
import dc_analysis, ticker, options, circuit, devices, printing, utilities, results

def ac_analysis(circ, start, nsteps, stop, step_type, xop=None, mna=None,\
	AC=None, Nac=None, J=None, data_filename="stdout", verbose=3):
	"""Performs an AC analysis of the circuit (described by circ).
	"""
	
	if data_filename == 'stdout':
		verbose = 0

	#check step/start/stop parameters
	if start == 0:
		printing.print_general_error("AC analysis has start frequency = 0")
                sys.exit(5)
	if start > stop:
		printing.print_general_error("AC analysis has start > stop")
		sys.exit(1)
	if nsteps < 1:
		printing.print_general_error("AC analysis has number of steps <= 1")
		sys.exit(1)
	if step_type == options.ac_log_step:
		omega_iter = utilities.log_axis_iterator(stop, start, nsteps)
	elif step_type == options.ac_lin_step:
		omega_iter = utilities.lin_axis_iterator(stop, start, nsteps)
	else:
		printing.print_general_error("Unknown sweep type.") 
		sys.exit(1)
	
	tmpstr = "Vea =", options.vea, "Ver =", options.ver, "Iea =", options.iea, "Ier =", \
	options.ier, "max_ac_nr_iter =", options.ac_max_nr_iter
	printing.print_info_line((tmpstr, 5), verbose)
	del tmpstr
	
	printing.print_info_line(("Starting AC analysis: ", 3), verbose)
	tmpstr = "w: start = %g Hz, stop = %g Hz, %d steps" % (start, stop, nsteps)
	printing.print_info_line((tmpstr, 3), verbose)
	del tmpstr

	#It's a good idea to call AC with prebuilt MNA matrix if the circuit is big
	if mna is None:
		(mna, N) = dc_analysis.generate_mna_and_N(circ)
		del N
		mna = utilities.remove_row_and_col(mna)
	if Nac is None:
		Nac = generate_Nac(circ)
		Nac = utilities.remove_row(Nac, rrow=0)
	if AC is None:
		AC = generate_AC(circ, [mna.shape[0], mna.shape[0]])
		AC = utilities.remove_row_and_col(AC)

	
	if circ.is_nonlinear():
		if J is not None:
			pass
			# we used the supplied linearization matrix
		else:
			if xop is None:
				printing.print_info_line(("Starting OP analysis to get a linearization point...", 3), verbose, print_nl=False)
				#silent OP
				xop = dc_analysis.op_analysis(circ, verbose=0)
				if xop is None: #still! Then op_analysis has failed!
					printing.print_info_line(("failed.", 3), verbose)
					printing.print_general_error("OP analysis failed, no linearization point available. Quitting.") 
					sys.exit(3)
				else:
					printing.print_info_line(("done.", 3), verbose)
			printing.print_info_line(("Linearization point (xop):", 5), verbose)
			if verbose > 4: xop.print_short()
			printing.print_info_line(("Linearizing the circuit...", 5), verbose, print_nl=False)
			J = generate_J(xop=xop.asmatrix(), circ=circ, mna=mna, Nac=Nac, data_filename=data_filename, verbose=verbose)
			printing.print_info_line((" done.", 5), verbose)
			# we have J, continue
	else: #not circ.is_nonlinear()
		# no J matrix is required.
		J = 0
	
	printing.print_info_line(("MNA (reduced):", 5), verbose)
	printing.print_info_line((str(mna), 5), verbose)
	printing.print_info_line(("AC (reduced):", 5), verbose)
	printing.print_info_line((str(AC), 5), verbose)
	printing.print_info_line(("J (reduced):", 5), verbose)
	printing.print_info_line((str(J), 5), verbose)
	printing.print_info_line(("Nac (reduced):", 5), verbose)
	printing.print_info_line((str(Nac), 5), verbose)
	
	sol = results.ac_solution(circ, ostart=start, ostop=stop, opoints=nsteps, stype=step_type, op=xop, outfile=data_filename)

	# setup the initial values to start the iteration:
	nv = len(circ.nodes_dict)
	j = numpy.complex('j')

	Gmin_matrix = dc_analysis.build_gmin_matrix(circ, options.gmin, mna.shape[0], verbose)

	iter_n = 0  # contatore d'iterazione
	#printing.print_results_header(circ, fdata, print_int_nodes=options.print_int_nodes, print_omega=True)
	printing.print_info_line(("Solving... ", 3), verbose, print_nl=False)
	tick = ticker.ticker(increments_for_step=1)
	tick.display(verbose > 1)

	x = xop
	for omega in omega_iter:
		(x, error, solved, n_iter) = dc_analysis.dc_solve(mna=(mna + numpy.multiply(j*omega, AC) + J), \
		Ndc=Nac,  Ntran=0, circ=circuit.circuit(title="Dummy circuit for AC", filename=None), Gmin=Gmin_matrix, x0=x, \
		time=None, locked_nodes=None, MAXIT=options.ac_max_nr_iter, skip_Tt=True, verbose=0)
		if solved:
			tick.step(verbose > 1)
			iter_n = iter_n + 1
			# hooray!
			sol.add_line(omega, x)
		else:
			break
	
	tick.hide(verbose > 1)
	
	if solved:
		printing.print_info_line(("done.", 3), verbose)
		ret_value = sol
	else:
		print "failed."
		ret_value =  None
	
	return ret_value

def generate_AC(circ, shape):
	"""Generates the AC coefficients matrix. 
	Shape is the REDUCED MNA shape, AC will be of the same shape.
	
	It's easy to set up the voltage lines, we know that line 2 refers to 
	node 2, etc... 
	
	A capacitor between nodes n1 and n2 determines the following elements:
	
	(KCL node n1) +j*w*C V(n1) - j*w*C V(n2) + ... = ...
	(KCL node n2) -j*w*C V(n1) + j*w*C V(n2) + ... = ...
	
	Inductors generate, together with voltage sources, ccvs, vcvs, a 
	additional line in the mna matrix, and hence in AC too. The current 
	flowing through the device gets added to the x vector.
	
	In inductors, we have:
	
	(KVL over n1 and n2) V(n1) - V(n2) - j*w*L I(inductor) = 0

	To understand on which line is the KVL line for an inductor, we use the 
	*order* of the elements in circuit.elements:
	First are all voltage lines, then the current ones in the same order of 
	the elements that introduce them.
	
	Returns: the UNREDUCED AC matrix
	"""
	AC = numpy.matrix(numpy.zeros((shape[0]+1, shape[1]+1)))
	nv = len(circ.nodes_dict)# - 1
	i_eq = 0 #each time we find a vsource or vcvs or ccvs, we'll add one to this.
	for elem in circ.elements:
		if isinstance(elem, devices.vsource) or isinstance(elem, devices.evsource) or \
		isinstance(elem, devices.hvsource):
			#notice that hvsources aren't yet implemented now!
			i_eq = i_eq + 1
		elif isinstance(elem, devices.capacitor):
			n1 = elem.n1
			n2 = elem.n2
			AC[n1, n1] = AC[n1, n1] + elem.C
			AC[n1, n2] = AC[n1, n2] - elem.C
			AC[n2, n2] = AC[n2, n2] + elem.C
			AC[n2, n1] = AC[n2, n1] - elem.C
		elif isinstance(elem, devices.inductor):
			AC[nv + i_eq, nv + i_eq] = -1 * elem.L
			if len(elem.coupling_devices):
				for cd in elem.coupling_devices:
					# get id+descr of the other inductor (eg. "L32")
					other_id_wdescr = cd.get_other_inductor("L"+elem.descr)
					# find its index to know which column corresponds to its current
					other_index = circ.find_vde_index(other_id_wdescr)
					# add the term.
					AC[nv + i_eq, nv + other_index] += -1 * cd.M
			i_eq = i_eq + 1
		
	if options.cmin > 0:
		cmin_mat = numpy.matrix(numpy.eye(shape[0]+1-i_eq))
		cmin_mat[0, 1:] = 1
		cmin_mat[1:, 0] = 1
		cmin_mat[0, 0] = cmin_mat.shape[0]-1
		AC[:-i_eq, :-i_eq] += options.cmin*cmin_mat

	return AC

def generate_Nac(circ):
	"""Generate the vector holding the contribution of AC sources.
	"""
	n_of_nodes = len(circ.nodes_dict)
	Nac = numpy.mat(numpy.zeros((n_of_nodes, 1)), dtype=complex)
	j = numpy.complex('j')
	# process isources
	for elem in circ.elements:
		if isinstance(elem, devices.isource) and elem.abs_ac is not None:
			#convenzione normale!
			N[elem.n1, 0] = N[elem.n1, 0] + elem.abs_ac*numpy.exp(j*elem.arg_ac)
			N[elem.n2, 0] = N[elem.n2, 0] - elem.abs_ac*numpy.exp(j*elem.arg_ac)
	# process vsources
	# for each vsource, introduce a new variable: the current flowing through it.
	# then we introduce a KVL equation to be able to solve the circuit
	for elem in circ.elements:
		if circuit.is_elem_voltage_defined(elem):
			index = Nac.shape[0] 
			Nac = utilities.expand_matrix(Nac, add_a_row=True, add_a_col=False)
			if isinstance(elem, devices.vsource) and elem.abs_ac is not None:
				Nac[index, 0] = -1.0*elem.abs_ac*numpy.exp(j*elem.arg_ac)
	return Nac

def generate_J(xop, circ, mna, Nac, data_filename, verbose=0):
	# setup J
	# build the linearized matrix (stored in J)
	J = numpy.mat(numpy.zeros(mna.shape))
	Tlin = numpy.mat(numpy.zeros(Nac.shape))
        for elem in circ.elements:
		if elem.is_nonlinear:
			dc_analysis.update_J_and_Tx(J, Tlin, xop, elem, time=None)
	#del Tlin # not needed! **DC**!

	return J

