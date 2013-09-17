# -*- coding: iso-8859-1 -*-
# shooting.py
# Shooting analysis module
# Copyright 2009 Giuseppe Venturini

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

"""Periodic steady state analysis based on the shooting method"""

import sys
import numpy, numpy.linalg
import transient, implicit_euler, dc_analysis, ticker, options, circuit, printing, utilities
import results, devices

def shooting(circ, period, step=None, x0=None, points=None, autonomous=False,
             mna=None, Tf=None, D=None, outfile='stdout', vector_norm=lambda v: max(abs(v)), verbose=3):
	"""Performs a periodic steady state analysis based on the algorithm described in
	Brambilla, A.; D'Amore, D., "Method for steady-state simulation of 
	strongly nonlinear circuits in the time domain," Circuits and 
	Systems I: Fundamental Theory and Applications, IEEE Transactions on, 
	vol.48, no.7, pp.885-889, Jul 2001
	URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=933329&isnumber=20194
	
	The results have been computed again by me, the formulas are not exactly the 
	same, but the idea behind the shooting algorithm is.
	
	This method allows us to have a period with many points without having to
	invert a huge matrix (and being limited to the maximum matrix size).

	A tran is performed to initialize the solver.
	
	We compute the change in the last point, calculating several matrices in
	the process.
	From that, with the same matrices we calculate the changes in all points, 
	starting from 0 (which is the same as the last one), then 1, ...

	Key points:
	- Only not autonomous circuits are supported.
	- The time step is constant
	- Implicit euler is used as DF
	
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
	autonomous has to be False, autonomous circuits are not supported
	outfile is the output filename. Defaults to stdout.
	verbose is set to zero (print errors only) if datafilename == 'stdout'.

	Returns: nothing
	"""
	
	if outfile == "stdout": 
		verbose = 0

	printing.print_info_line(("Starting periodic steady state analysis:",3), verbose)
	printing.print_info_line(("Method: shooting",3), verbose)
	
	if isinstance(x0, results.op_solution):
		x0 = x0.asmatrix()
	if mna is None or Tf is None:
		(mna, Tf) = dc_analysis.generate_mna_and_N(circ)
		mna = utilities.remove_row_and_col(mna)
		Tf = utilities.remove_row(Tf, rrow=0)
	elif not mna.shape[0] == Tf.shape[0]:
		printing.print_general_error("mna matrix and N vector have different number of rows.")
		sys.exit(0)
	
	if D is None:
		D = transient.generate_D(circ, [mna.shape[0], mna.shape[0]])
		D = utilities.remove_row_and_col(D)
	elif not mna.shape == D.shape:
		printing.print_general_error("mna matrix and D matrix have different sizes.")
		sys.exit(0)
	
	(points, step) = check_step_and_points(step, points, period)
	print "points", points
	print "step", step
	
	n_of_var = mna.shape[0]
	locked_nodes = circ.get_locked_nodes()
	
	printing.print_info_line(("Starting transient analysis for algorithm init: tstop=%g, tstep=%g... " % (10*points*step, step),3), verbose, print_nl=False)
	xtran = transient.transient_analysis(circ=circ, tstart=0, tstep=step, tstop=10*points*step, method="TRAP", x0=None, mna=mna, N=Tf, \
        D=D, use_step_control=False, outfile=outfile+".tran", return_req_dict={"points":points}, verbose=0)
	if xtran is None:
		print "failed."
		return None
	printing.print_info_line(("done.",3), verbose)
	
	x = []
	for index in range(points):
		x.append(xtran[index*n_of_var:(index+1)*n_of_var,0])

	tick = ticker.ticker(increments_for_step=1)

	MAass_static, MBass = build_static_MAass_and_MBass(mna, D, step)
	
	# This contains
	# the time invariant part, Tf
	# time variable component: Tt this is always the same, since the time interval is the same
	# this holds all time-dependent sources (both V/I).
	Tass_static_vector = build_Tass_static_vector(circ, Tf, points, step, tick, n_of_var, verbose)
	
	converged = False
	printing.print_info_line(("Solving... ",3), verbose, print_nl=False)
	tick.reset()
	tick.display(verbose > 2)
	
	iteration = 0 # newton iteration counter
	conv_counter = 0

	while True:
		dx = []
		Tass_variable_vector = []
		MAass_variable_vector = []
		for index in range(points):
			if index == 0:
				xn_minus_1 = x[points-1]
			else:
				xn_minus_1 = x[index-1]
			MAass_variable, Tass_variable = get_variable_MAass_and_Tass(circ, x[index], xn_minus_1, mna, D, step, n_of_var)
			MAass_variable_vector.append(MAass_variable + MAass_static)
			Tass_variable_vector.append(Tass_variable + Tass_static_vector[index])
		
		dxN = compute_dxN(circ, MAass_variable_vector, MBass, Tass_variable_vector, n_of_var, points, verbose=verbose)
		td = dc_analysis.get_td(dxN, locked_nodes, n=-1)
		x[points-1] = td * dxN + x[points-1]

		for index in range(points-1):
			if index == 0:
				dxi_minus_1 = dxN
			else:
				dxi_minus_1 = dx[index-1]
			dx.append(compute_dx(MAass_variable_vector[index], MBass, Tass_variable_vector[index], dxi_minus_1))
			td = dc_analysis.get_td(dx[index], locked_nodes, n=-1)
			x[index] = td*dx[index] + x[index]
		dx.append(dxN)

		if (vector_norm_wrapper(dx, vector_norm) < min(options.ver, options.ier)*vector_norm_wrapper(x,vector_norm) + min(options.vea, options.iea)): #\
		#and (dc_analysis.vector_norm(residuo) < options.er*dc_analysis.vector_norm(x) + options.ea):
			if conv_counter ==3:
				converged = True
				break
			else:
				conv_counter = conv_counter + 1
		elif vector_norm(dx[points-1]) is numpy.nan: #needs work fixme
			raise OverflowError
			#break
		else:
			conv_counter = 0
			tick.step(verbose > 2)
	
		if options.shooting_max_nr_iter and iteration == options.shooting_max_nr_iter:
			printing.print_general_error("Hitted SHOOTING_MAX_NR_ITER (" + str(options.shooting_max_nr_iter) + "), iteration halted.")
			converged = False
			break
		else:
			iteration = iteration + 1

	tick.hide(verbose > 2)
	if converged:
		printing.print_info_line(("done.", 3), verbose)
		t = numpy.mat(numpy.arange(points)*step)
		t = t.reshape((1, points))
		xmat = x[0]
		for index in xrange(1, points):
			xmat = numpy.concatenate((xmat, x[index]), axis=1)
		sol = results.pss_solution(circ=circ, method="shooting", period=period, outfile=outfile, t_array=t, x_array=xmat)
		#print_results(circ, x, fdata, points, step)
	else:
		print "failed."
		sol = None
	return sol

def vector_norm_wrapper(vector, norm_fun):
	max = 0
	for elem in vector:
		new_max = norm_fun(elem)
		if max < new_max:
			max = new_max
	return max


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
	
	return (points, step)

def build_static_MAass_and_MBass(mna, D, step):
	(C1, C0) = implicit_euler.get_df_coeff(step)
	MAass = mna + D*C1
        MBass = D * C0
	return (MAass, MBass)

def build_Tass_static_vector(circ, Tf, points, step, tick, n_of_var, verbose=3):
        Tass_vector = []
	nv = len(circ.nodes_dict)
	printing.print_info_line(("Building Tass...", 5), verbose, print_nl=False)

        tick.reset()
        tick.display(verbose > 2)
        for index in xrange(0, points):
                Tt = numpy.zeros((n_of_var, 1))
                v_eq = 0
                time = index * step
                for elem in circ.elements:
                        if (isinstance(elem, devices.vsource) or isinstance(elem, devices.isource)) and elem.is_timedependent:
                                if isinstance(elem, devices.vsource):
                                        Tt[nv - 1 + v_eq, 0] = -1.0 * elem.V(time)
                                elif isinstance(elem, devices.isource):
                                        if elem.n1:
                                                Tt[elem.n1-1, 0] = \
                                                Tt[elem.n1-1, 0] + elem.I(time)
                                        if elem.n2:
                                                Tt[elem.n2-1, 0] = \
                                                Tt[elem.n2-1, 0] - elem.I(time)
                        if circuit.is_elem_voltage_defined(elem):
                                v_eq = v_eq +1
                tick.step(verbose > 2)
                Tass_vector.append(Tf+Tt)
        tick.hide(verbose > 2)
	printing.print_info_line(("done.", 5), verbose)

        return Tass_vector

def get_variable_MAass_and_Tass(circ, xi, xi_minus_1, M, D, step, n_of_var):
	Tass = numpy.zeros((n_of_var, 1))
	J = numpy.zeros((n_of_var, n_of_var))
	(C1, C0) = implicit_euler.get_df_coeff(step)

	for elem in circ.elements:
		# build all dT(xn)/dxn (stored in J) and T(x)
		if elem.is_nonlinear:
			output_ports = elem.get_output_ports()
			for index in range(len(output_ports)):
				n1, n2 = output_ports[index]
				ports = elem.get_drive_ports(index)
				v_ports = []
				for port in ports:
					v = 0 # build v: remember we trashed the 0 row and 0 col of mna -> -1
					if port[0]:
						v = v + xi[port[0] - 1, 0]
					if port[1]:
						v = v - xi[port[1] - 1, 0]
					v_ports.append(v)
				if n1:
					Tass[n1 - 1, 0] = Tass[n1 - 1, 0] + elem.i(index, v_ports)
				if n2:
					Tass[n2 - 1, 0] = Tass[n2 - 1, 0] - elem.i(index, v_ports)
				for pindex in xrange(len(ports)):
					if n1:
						if ports[pindex][0]:
							J[n1 - 1, ports[pindex][0]-1] = \
							J[n1 - 1, ports[pindex][0]-1] + elem.g(index, v_ports, pindex)
						if ports[pindex][1]:
							J[n1 - 1, ports[pindex][1]-1] =\
							J[n1 - 1, ports[pindex][1]-1] - 1.0*elem.g(index, v_ports, pindex)
					if n2:
						if ports[pindex][0]:
							J[n2 - 1, ports[pindex][0]-1] = \
							J[n2 - 1, ports[pindex][0]-1] - 1.0*elem.g(index, v_ports, pindex)
						if ports[pindex][1]:
							J[n2 - 1, ports[pindex][1]-1] =\
							J[n2 - 1, ports[pindex][1]-1] + elem.g(index, v_ports, pindex)
	
	Tass = Tass + D*C1*xi + M*xi + D*C0*xi_minus_1

	return (J, Tass)

def compute_dxN(circ, MAass_vector, MBass, Tass_vector, n_of_var, points, verbose=3):
	temp_mat1 = numpy.mat(numpy.eye(n_of_var))
	for index in range(points):
		temp_mat1 = -1*numpy.linalg.inv(MAass_vector[index])*MBass*temp_mat1
	temp_mat2 = numpy.mat(numpy.zeros((n_of_var,1)))
	for index in range(points):
		temp_mat3 = -1*numpy.linalg.inv(MAass_vector[index])*Tass_vector[index]
		for index2 in range(index+1, points):
			temp_mat3 = -1*numpy.linalg.inv(MAass_vector[index2])*MBass*temp_mat3
		temp_mat2 = temp_mat2 + temp_mat3

	dxN = numpy.linalg.inv(numpy.mat(numpy.eye(n_of_var)) - temp_mat1) * temp_mat2

	return dxN

def compute_dx(MAass, MBass, Tass, dxi_minus_1):
	dxi = -1 * numpy.linalg.inv(MAass) * (MBass * dxi_minus_1 + Tass)
	return dxi

