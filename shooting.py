# -*- coding: iso-8859-1 -*-
# shooting.py
# Shooting analysis module
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

"""Shooting analysis module"""

import sys
import numpy
import transient, implicit_euler, dc_analysis, ticker, options, circuit, printing, utilities

#_debug = True

def shooting(circ, period, step=None, mna=None, Tf=None, D=None, points=None, autonomous=False, x0=None,  data_filename='stdout', vector_norm=lambda v: max(abs(v)), verbose=3):
	"""Performs a shooting analysis. 
	
	Time step is constant, IE will be used as DF
	
	Parameters:
	circ is the circuit description class
	period is the period of the solution
	mna, D, Tf are not compulsory they will be computed if they're set to None
	step is the time step between consecutive points
	points is the number of points to be used
	step and points are mutually exclusive options:
	- if step is specified, the number of points will be automatically determined 
	- if points is set, the step will be automatically determined
	- if none of them is set, options.shooting_default_points will be used as points
	autonomous has to be False, not autonomous circuits are not supported
	x0 is the initial guess to be used. Needs work.
	datafilename is the output filename. Defaults to stdout.
	
	Returns: nothing
	"""
	if data_filename != "stdout" and verbose > 2:
		print "Starting shooting analysis:"
	
	if data_filename != "stdout":
		fdata = open(data_filename, "w")
	else:
		fdata = sys.stdout
	
	if mna is None or Tf is None:
		(mna, Tf) = dc_analysis.generate_mna_and_N(circ)
		mna = utilities.remove_row_and_col(mna)
		Tf = utilities.remove_row(Tf, rrow=0)
	elif not mna.shape[1] == Tf.shape[1]:
		printing.print_general_error("mna matrix and N vector have different number of columns.")
		sys.exit(0)
	
	if D is None:
		D = transient.generate_D(circ, [mna.shape[0], mna.shape[0]])
		D = utilities.remove_row_and_col(D)
	elif not mna.shape == D.shape:
		printing.print_general_error("mna matrix and D matrix have different sizes.")
		sys.exit(0)
	
	(points, step) = check_step_and_points(step, points, period)
	
	#nv = len(circ.nodes_dict)	
	n_of_var = mna.shape[0]
	locked_nodes = circ.get_locked_nodes()
	if verbose > 2: tick = ticker.ticker(increments_for_step=1)

	CMAT = build_CMAT(mna, D, step, points, tick, n_of_var=n_of_var, \
		print_out=(data_filename!="stdout"), verbose=verbose)

	x = build_x(mna, step, points, tick, x0=x0, n_of_var=n_of_var, \
		print_out=(data_filename!="stdout"), verbose=verbose)

	Tf = build_Tf(Tf, points, tick, n_of_var=n_of_var, print_out=(data_filename!="stdout"), verbose=verbose)
	
	# time variable component: Tt this is always the same in each iter. So we build it once for all
	# this holds all time-dependent sources (both V/I).
	Tt = build_Tt(circ, points, step, tick, n_of_var=n_of_var, print_out=(data_filename!="stdout"), verbose=verbose)
	
	converged = False
	if verbose > 2: 
		if data_filename != "stdout": 
			sys.stdout.write("Solving... ")
		tick.reset()
		tick.display()
	J = numpy.mat(numpy.zeros(CMAT.shape))
	T = numpy.mat(numpy.zeros((CMAT.shape[0], 1)))
	td = numpy.mat(numpy.zeros((points, 1)))
	iteration = 0 # newton iteration counter
	
	while True:
		if iteration: # the first time are already all zeros
			J[:, :] = 0
			T[:, 0] = 0
			td[:, 0] = 0
		for index in xrange(points):
			for elem in circ.elements:
				# build all dT(xn)/dxn (stored in J) and T(x)
				if elem.is_nonlinear:
					ports = elem.get_ports()
					v_ports = []
					for port in ports:
						v = 0 # build v: remember we trashed the 0 row and 0 col of mna -> -1
						if port[0]:
							v = v + x[index*n_of_var + port[0] - 1, 0]
						if port[1]:
							v = v - x[index*n_of_var + port[1] - 1, 0]
						v_ports.append(v)
					if elem.n1:
						T[index*n_of_var + elem.n1 - 1, 0] = T[index*n_of_var + elem.n1 - 1, 0] + elem.i(v_ports)
					if elem.n2:
						T[index*n_of_var + elem.n2 - 1, 0] = T[index*n_of_var + elem.n2 - 1, 0] - elem.i(v_ports)
					for pindex in xrange(len(ports)):
						if elem.n1:
							if ports[pindex][0]:
								J[index*n_of_var + elem.n1-1, index*n_of_var + ports[pindex][0]-1] = \
								J[index*n_of_var + elem.n1-1, index*n_of_var + ports[pindex][0]-1] + elem.g(v_ports, pindex)
							if ports[pindex][1]:
								J[index*n_of_var + elem.n1-1, index*n_of_var + ports[pindex][1]-1] =\
								J[index*n_of_var + elem.n1-1, index * n_of_var + ports[pindex][1]-1] - 1.0*elem.g(v_ports, pindex)
						if elem.n2:
							if ports[pindex][0]:	
								J[index*n_of_var + elem.n2-1, index*n_of_var + ports[pindex][0]-1] = \
								J[index*n_of_var + elem.n2-1, index*n_of_var + ports[pindex][0]-1] - 1.0*elem.g(v_ports, pindex)
							if ports[pindex][1]:
								J[index*n_of_var + elem.n2-1, index*n_of_var + ports[pindex][1]-1] =\
								J[index*n_of_var + elem.n2-1, index*n_of_var + ports[pindex][1]-1] + elem.g(v_ports, pindex)

		J = J + CMAT
		residuo = CMAT*x + T + Tf + Tt
		dx = numpy.linalg.inv(J) * (-1 * residuo)
		#td
		for index in xrange(points):
			td[index, 0] = dc_analysis.get_td(dx[index*n_of_var:(index+1)*n_of_var, 0], locked_nodes, n=-1)
		x = x + min(abs(td))[0, 0] * dx
		#print x[points*n_of_var, 0]
		# convergence and maxit test FIXME: WHAT?
		if (vector_norm(dx) < min(options.ver, options.ier)*vector_norm(x) + min(options.vea, options.iea)): #\
		#and (dc_analysis.vector_norm(residuo) < options.er*dc_analysis.vector_norm(x) + options.ea):
			converged = True
			break
		elif vector_norm(dx) is numpy.nan: #needs work fixme
			raise OverflowError
			#break
		else:
			if verbose > 2: 
				tick.step()
	
		if options.shooting_max_nr_iter and iteration == options.shooting_max_nr_iter:
			printing.print_general_error("Hitted SHOOTING_MAX_NR_ITER (" + str(options.shooting_max_nr_iter) + "), iteration halted.")
			converged = False
			break
		else:
			iteration = iteration + 1
	if verbose > 2: 
		tick.hide()
	if converged:
		if verbose > 2 and data_filename != "stdout": 
			print "done."
		printing.print_results_header(circ, fdata, print_int_nodes=options.print_int_nodes, print_time=True)
		

		for index in xrange(points):
			printing.print_results_on_a_line(time=index*step, x=x[index*n_of_var:(index+1)*n_of_var, 0], fdata=fdata, circ=circ, print_int_nodes=options.print_int_nodes, iter_n=0)
			#printing.print_results_at_time(index*step, x[index*n_of_var:(index+1)*n_of_var, 0], fdata, iter_n=0)
	else:
		if verbose > 2 and data_filename != "stdout": 
			print "failed."
	return


def set_submatrix(row, col, dest_matrix, source_matrix):
	"""Copies the source_matrix in dest_matrix, 
	the position of the upper left corner of source matrix is (row, col) within dest_matrix
	
	Returns dest_matrix
	"""
	for li in xrange(source_matrix.shape[0]):
		for ci in xrange(source_matrix.shape[1]):
			if source_matrix[li, ci] != 0:
				dest_matrix[row + li, col + ci] = source_matrix[li, ci]
	return dest_matrix
	
def get_e(index, length):
	"""Builds a e(j=index) col vector
	e(j) is defined as:
		e(i, 0) = 1 if i=j
		e(i, 0) = 0 otherwise
	
	Returns: e(index)
	"""
	e = numpy.mat(numpy.zeros((length, 1)))
	e[index, 0] = 1
	return e
	
def check_step_and_points(step, points, period):
	"""Sets consistently the step size and the number of points, according to the given period
	Returns: (points, step)
	"""
	if step is None and points is None:
		print "Warning: shooting had no step nor n. of points setted. Using", options.shooting_default_points,"points."
		points = options.shooting_default_points
	elif step is not None and points is not None:
		print "Warning: shooting had both step and n. of points setted. Using", step, "step. (NA)"
		points = None
		
	if points:
		step = (1.0 * period)/(points - 1)
	else:
		points = (1.0 * period) / step
		if points % 1 != 0:
			step = step + (step * (points % 1)) / int(points)
			points = int((1.0 * period) / step)
			print "Warning: adapted step is", step
		else:
			points = int(points)
		points = points + 1 #0 - N where xN is in reality the first point of the second period!!
	
	return (points, step)

def build_CMAT(mna, D, step, points, tick, n_of_var=None, print_out=True, verbose=3):
	if n_of_var is None: n_of_var = mna.shape[0]
	if verbose > 4 and print_out: 
		sys.stdout.write("Building the CMAT ("+str(n_of_var*points)+"x"+str(n_of_var*points)+")... ")
	if verbose > 2: 
		tick.reset()
		tick.display()
	(C1, C0) = implicit_euler.get_df_coeff(step)
	I = numpy.mat(numpy.eye(n_of_var))
	M = mna + C1*D
	N = C0 * D
	#Z = numpy.mat(numpy.zeros((n_of_var, n_of_var)))
	CMAT = numpy.mat(numpy.zeros((n_of_var*points, n_of_var*points)))
	for li in xrange(points): #li = line index
		for ci in xrange(points):
			if li == 0:
				if ci == 0:
					temp = 1.0 * I
				elif ci == points - 1:
					temp = -1.0 * I
				else:
					continue #temp = Z
			else:
				if ci == li:
					temp = M
				elif ci == li -1:
					temp = N
				else:
					continue #temp = Z
			CMAT = set_submatrix(row=li*n_of_var, col=ci*n_of_var, dest_matrix=CMAT, source_matrix=temp)
		if verbose > 2: tick.step()
	if verbose > 2: tick.hide()
	if verbose > 4 and print_out: print "done."
	#print CMAT
	return CMAT
	
def build_x(mna, step, points, tick, x0=None, n_of_var=None, print_out=True, verbose=3):
	if n_of_var is None: n_of_var = mna.shape[0]
	if verbose > 4 and print_out: sys.stdout.write("Building x...")
	if verbose > 2: 
		tick.reset()
		tick.display()
	x = numpy.mat(numpy.zeros((points*n_of_var, 1)))
	if x0 is not None:
		if x0.shape[0] != n_of_var:
			print "Warning x0 has the wrong dimensions. Using all-0."
		else:
			for index in xrange(points):
				x = set_submatrix(row=index*n_of_var, col=0, dest_matrix=x, source_matrix=x0)
				if verbose > 2: tick.step()
	if verbose > 2: tick.hide()
	if verbose > 4 and print_out: print "done."
	
	return x
	
def build_Tf(sTf, points, tick, n_of_var, print_out=True, verbose=3):
	if verbose > 4 and print_out: sys.stdout.write("Building Tf... ")
	if verbose > 2: 
		tick.reset()
		tick.display()
	Tf = numpy.mat(numpy.zeros((points*n_of_var, 1)))
	for index in xrange(points):
		Tf = set_submatrix(row=index*n_of_var, col=0, dest_matrix=Tf, source_matrix=sTf)
		if verbose > 2: tick.step()
	if verbose > 2: tick.hide()
	if verbose > 4 and print_out: print "done."
	
	return Tf
	
def build_Tt(circ, points, step, tick, n_of_var, print_out=True, verbose=3):
	nv = len(circ.nodes_dict)
	if verbose > 4 and print_out: sys.stdout.write("Building Tt... ")
	if verbose > 2: 
		tick.reset()
		tick.display()	
	Tt = numpy.zeros((points*n_of_var, 1))
	for index in xrange(points):
		v_eq = 0
		time = index * step
		for elem in circ.elements:
			if (isinstance(elem, circuit.vsource) or isinstance(elem, circuit.isource)) and elem.is_timedependent:
				if isinstance(elem, circuit.vsource):
					Tt[index*n_of_var + nv - 1 + v_eq, 0] = -1.0 * elem.V(time)
				elif isinstance(elem, circuit.isource):
					if elem.n1:
						Tt[index*n_of_var + elem.n1-1, 0] = \
						Tt[index*n_of_var + elem.n1-1, 0] + elem.I(time)
					if elem.n2:
						Tt[index*n_of_var + elem.n2-1, 0] = \
						Tt[index*n_of_var + elem.n2-1, 0] - elem.I(time)
			if circuit.is_elem_voltage_defined(elem):
				v_eq = v_eq +1
			#print Tt[index*n_of_var:(index+1)*n_of_var]
		if verbose > 2: tick.step()
	if verbose > 2: tick.hide()
	if verbose > 4 and print_out: print "done."
	
	return Tt
