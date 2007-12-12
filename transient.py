# -*- coding: iso-8859-1 -*-
# transient.py
# Transient analysis
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

""" This module offers the methods required to perform a transient analysis.
Our problem can be written as:
	D*dx/dt + MNA*x + Tv(x) + Tt(t) + N = 0
We need:
	1. the mna matrix MNA
	2. N
	3. T(x)
	4. Tt(t) has to be evaluated step by step
	5. D matrix
	6. a differentiation method to approximate dx/dt
"""

import sys, imp
import numpy
import dc_analysis, implicit_euler, ticker, options, circuit, printing, utilities


#methods, add here
IMPLICIT_EULER = "IMPLICIT_EULER"
TRAP = "TRAP"
GEAR1 = "GEAR1"
GEAR2 = "GEAR2"
GEAR3 = "GEAR3"
GEAR4 = "GEAR4"
GEAR5 = "GEAR5"
GEAR6 = "GEAR6"

def transient_analysis(circ, tstart, tstep, tstop, method=TRAP, x0=None, mna=None, N=None, \
	data_filename="stdout", use_step_control=True,  \
	print_step_and_lte=False, verbose=3):
	"""Performs a transient analysis of the circuit described by circ.
	
	Important parameters:
	- tstep is the maximum step to be allowed during simulation.
	- print_step_and_lte is a boolean value. Is set to true, the step and the LTE of the first
	element of x will be printed out to step_and_lte.graph in the current directory.
	
	"""
	#use_aposteriori_step_control = False
	#use_step_control = False
	#print_step_and_lte = True
	
	HMAX = tstep
	
	#check paramenters
	if tstart > tstop:
		printing.print_general_error("tstart > tstop")
		sys.exit(1)
	if tstep < 0:
		printing.print_general_error("tstep < 0")
		sys.exit(1)

	if verbose > 4:
		print "Vea =", options.vea, "Ver =", options.ver, "Iea =", options.iea, "Ier =", \
		options.ier, "max_time_iter =", options.transient_max_time_iter, "HMIN=",options.hmin
	#	print "tstart = "+str(tstart), "tstop = " + str(tstop), "tstep = "+str(tstep)
	
	locked_nodes = circ.get_locked_nodes()
	
	# setup output streams:
	if data_filename != "stdout":
		fdata = open(data_filename, "w")
	else:
		fdata = sys.stdout
	if print_step_and_lte:
		flte = open("step_and_lte.graph", "w")
		flte.write("#T\tStep\tLTE\n")
	
	if verbose > 2 and fdata is not sys.stdout:
		print "Starting transient analysis: "
		print "Selected method: "+method
	#It's a good idea to call transient with prebuilt MNA and N matrix
	if mna is None or N is None:
		(mna, N) = dc_analysis.generate_mna_and_N(circ)
		mna = utilities.remove_row_and_col(mna)
		N = utilities.remove_row(N, rrow=0)
	elif not mna.shape[1] == N.shape[1]:
		printing.print_general_error("mna matrix and N vector have different number of columns.")
		sys.exit(0)
	
	# We generate the D matrix by default. It could be provided but that makes sense only 
	# if you do more than one tran analysis. 
	# Ahkab is not ready to that. For example, output streams should be changed...
	D = generate_D(circ, [mna.shape[0], mna.shape[0]])
	D = utilities.remove_row_and_col(D)

	# setup x0
	if x0 is None:
		if verbose > 4 and fdata is not sys.stdout: 
			print "Generating x(t="+str(tstart)+") = 0"
		x0 = numpy.matrix(numpy.zeros((mna.shape[0], 1)))
	else:
		if verbose > 4 and fdata is not sys.stdout:
			print "Using the supplied op as x(t="+str(tstart)+")."

	if verbose > 4 and fdata is not sys.stdout:
		print "x0:"
		printing.print_results_header(circ, sys.stdout, print_int_nodes=True, print_time=False)
		printing.print_results_on_a_line(None, x0, sys.stdout, circ, print_int_nodes=True, iter_n=10)
		#print x0
	
	# setup the df method
	if verbose > 4: 
		sys.stdout.write("Selecting the appropriate DF ("+method+")... ")
	if method == IMPLICIT_EULER:
		import implicit_euler as df
	elif method == TRAP:
		import trap as df
	elif method == GEAR1:
		import gear as df
		df.order = 1
	elif method == GEAR2:
		import gear as df
		df.order = 2
	elif method == GEAR3:
		import gear as df
		df.order = 3
	elif method == GEAR4:
		import gear as df
		df.order = 4
	elif method == GEAR5:
		import gear as df
		df.order = 5
	elif method == GEAR6:
		import gear as df
		df.order = 6
	else:
		df = import_custom_df_module(method, print_out=(fdata != 'stdout'))
		# df is none if module is not found
	
	if df is None:
		sys.exit(23)
		
	if not df.has_ff() and use_step_control:
		printing.print_warning("The chosen DF does not support step control. Turning off the feature.")
		use_step_control = False
		#use_aposteriori_step_control = False
	elif verbose > 4:
		sys.stdout.write("done\n")
	
	# setup the data buffer
	# if you use the step control, the buffer has to be one point longer.
	# That's because the excess point is used by a FF in the df module to predict the next value.
	if verbose > 4:
		sys.stdout.write("Setting up the buffer... ")
	((max_x, max_dx), (pmax_x, pmax_dx)) = df.get_required_values()
	if max_x is None and max_dx is None:
		printing.print_general_error("df doesn't need any value?")
		sys.exit(1)
	if use_step_control:
		thebuffer = dfbuffer(length=max(max_x, max_dx, pmax_x, pmax_dx) + 1, width=3)
	else:
		thebuffer = dfbuffer(length=max(max_x, max_dx) + 1, width=3)
	thebuffer.add((tstart, x0, None)) #setup the fist values
	if verbose > 4: 
		sys.stdout.write("done\n")
	
	# import implicit_euler to be used in the first iterations
	# this is because we don't have any dx when we start, nor any past point value
	if (max_x is not None and max_x > 0) or max_dx is not None:
		import implicit_euler
	
	if verbose > 4:
		print "MNA (reduced):"
		print mna
		print "D (reduced):"
		print D
	
	# setup the initial values to start the iteration:
	x = None
	time = tstart
	nv = len(circ.nodes_dict)
	# lo step viene generato automaticamente, ma non superare mai quello fornito.
	if use_step_control:
		#tstep = min((tstop-tstart)/9999.0, HMAX, 100.0 * options.hmin)
		tstep = min((tstop-tstart)/9999.0, HMAX)
		if verbose > 4:
			print "Initial step:", tstep
	#else:
		#tstep = HMAX #should already be so, but, ynk
	if max_dx is None:
		max_dx_plus_1 = None
	else:
		max_dx_plus_1 = max_dx +1
	if pmax_dx is None:
		pmax_dx_plus_1 = None
	else:
		pmax_dx_plus_1 = pmax_dx +1
	
	# setup error vectors
	aerror = numpy.mat(numpy.zeros((x0.shape[0], 1)))
	aerror[:nv-1, 0] = options.vea
	aerror[nv-1:, 0] = options.vea
	rerror = numpy.mat(numpy.zeros((x0.shape[0], 1)))
	rerror[:nv-1, 0] = options.ver
	rerror[nv-1:, 0] = options.ier
	
	iter_n = 0  # contatore d'iterazione
	lte = None
	printing.print_results_header(circ, fdata, print_int_nodes=options.print_int_nodes, print_time=True)
	printing.print_results_on_a_line(time, x0, fdata, circ, print_int_nodes=options.print_int_nodes, iter_n=0)
	#printing.print_results_at_time(time, x0, fdata, iter_n)
	if fdata != sys.stdout:
		if verbose > 2:
			sys.stdout.write("Solving... ")
		if verbose > 1:
			tick = ticker.ticker(increments_for_step=1)
			tick.display()
	while time < tstop:
		if iter_n < max(max_x, max_dx_plus_1):
			[x_coeff, const, x_lte_coeff, prediction, pred_lte_coeff] = \
			implicit_euler.get_df([thebuffer.get_df_vector()[0]], tstep, \
			predict=(use_step_control and (iter_n >= max(pmax_x, pmax_dx_plus_1))))
		else:
			[x_coeff, const, x_lte_coeff, prediction, pred_lte_coeff] = \
			df.get_df(thebuffer.get_df_vector(), tstep, predict=use_step_control)
		
		if options.transient_prediction_as_x0 and use_step_control and prediction is not None:
			x0 = prediction
		elif x is not None:
			x0 = x
		
		(x1, error, solved) = dc_analysis.dc_solve(mna=(mna + numpy.multiply(x_coeff, D)) , N=(N + D*const), circ=circ, use_gmin=True, x0=x0, time=(time + tstep), locked_nodes=locked_nodes, MAXIT=options.transient_max_nr_iter, verbose=0)
		
		if solved:
			old_step = tstep #we will modify it, if we're using step control otherwise it's the same
			# step control (yeah)
			if use_step_control:
				if x_lte_coeff is not None and pred_lte_coeff is not None and prediction is not None:
					# this is the Local Truncation Error :)
					lte = abs((x_lte_coeff / (pred_lte_coeff - x_lte_coeff)) * (prediction - x1))
					# it should NEVER happen that new_step > 2*tstep, for stability
					new_step_coeff = 2 
					for index in xrange(x.shape[0]):
						if lte[index, 0] != 0:
							new_value = ((aerror[index, 0] + rerror[index, 0]*abs(x[index, 0])) / lte[index, 0]) \
							** (1.0 / (df.order+1))
							if new_value < new_step_coeff:
								new_step_coeff = new_value
							#print new_value
					new_step = tstep * new_step_coeff
					if options.transient_use_aposteriori_step_control and new_step < options.transient_aposteriori_step_threshold * tstep: 
						#don't recalculate a x for a small change
						tstep = check_step(new_step, time, tstop, HMAX)
						#print "Apost. (reducing) step = "+str(tstep)
						continue
					tstep = check_step(new_step, time, tstop, HMAX) # used in the next iteration
					#print "Apriori tstep = "+str(tstep)
				else:
					#print "LTE not calculated."
					lte = None
			if print_step_and_lte and lte is not None: 
				#if you wish to look at the step. We print just a lte
				flte.write(str(time)+"\t"+str(old_step)+"\t"+str(lte.max())+"\n")
			# if we get here, or aposteriori_step_control is disabled, or it's enabled and the error is small
			# enough. Anyway, the result is GOOD, STORE IT.
			time = time + old_step
			x = x1
			iter_n = iter_n + 1
			printing.print_results_on_a_line(time, x, fdata, circ, options.print_int_nodes, iter_n)
			#printing.print_results_at_time(time, x, fdata, iter_n)
			#thebuffer.add((numpy.mat(numpy.array(time)), x, numpy.multiply(x_coeff, x) + const))
			thebuffer.add((time, x, numpy.multiply(x_coeff, x) + const))
			if fdata != sys.stdout:
				if verbose > 1:
					tick.step()
		else:
			# If we get here, Newton failed to converge. We need to reduce the step...
			if use_step_control:
				tstep = tstep/5
				tstep = check_step(tstep, time, tstop, HMAX)
				if verbose > 4 and fdata is not sys.stdout:
					print "At "+str(time)+" reducing step: "+str(tstep)+" (convergence failed)"
			else: #we can't reduce the step
				printing.print_general_error("Can't converge with step "+str(tstep)+".")
				printing.print_general_error("Try setting --t-max-nr to a higher value or set step to a lower one.")
				solved = False
				break
		if options.transient_max_time_iter and iter_n == options.transient_max_time_iter:
			printing.print_general_error("Hitted MAX_TIME_ITER ("+str(options.transient_max_time_iter)+"), iteration halted.")
			solved = False
			break
	#end of while
	if not fdata == sys.stdout:
		fdata.close()
	
	if print_step_and_lte:
		flte.close()
	
	if fdata != sys.stdout and verbose > 1:
		tick.hide()
	
	if solved:
		if fdata != sys.stdout and verbose > 2:
			print "done."
			print "Average time step:", (tstop - tstart)/iter_n
	else:
		if fdata != sys.stdout:
			print "failed."
		(time, x) =  (None, None)
	
	return (time, x)

def check_step(tstep, time, tstop, HMAX):
	"""Checks the step for the following problems:
	- the step must be shorter than HMAX (that usually is the tstep provided by the user)
	- the step must be shorter than the simulation time left (ie tstop - time)
	- the step must be longer than options.hmin, if not halt the simulation.
	
	Returns: the step provided if it's ok, a shortened step otherwise.
	"""
	if tstep > HMAX:
		tstep = HMAX
	if tstop - time < tstep:
		tstep = tstop - time
	elif tstep < options.hmin:
		printing.print_general_error("Step size too small: "+str(tstep))
		raise Exception, "Step size too small"
	return tstep

def generate_D(circ, shape):
	"""Generates the derivate coefficients. Shape is the REDUCED MNA shape, D will be of the same shape.
	It's easy to set up the voltage lines, we know that line 2 refers to node 2, etc... 
	So everything's fine with capacitors. 
	Inductors generate, together with voltage sources, ccvs, vcvs, a additional line in the
	mna matrix, and hence in D too. The current flowing through the device gets added to the x vector.
	In inductors, we have:
	 V(n1) - V(n2) - VL = 0
	Where VL = L dI/dt
	That's 0 (zero) in DC analysis, but not in transient analysis, where it needs to be differentiated.
	To understand on which line does the inductor's L*dI/dt go, we use the order in circuit.elements:
	First are all voltage lines, then the current ones in the same order of the elements that introduce
	them.
	Therefore, we look at circ.elements.
	
	For every time t, the D matrix is used (elsewhere) to solve the following system:
	
	D*dx/dt + MNA*x + N  + T(x) = 0
	
	Returns: the UNREDUCED D matrix
	"""
	shape[0] = shape[0] + 1
	shape[1] = shape[1] + 1
	D = numpy.matrix(numpy.zeros(shape))
	nv = len(circ.nodes_dict)# - 1
	i_eq = 0 #each time we find a vsource or vcvs or ccvs, we'll add one to this.
	for elem in circ.elements:
		if isinstance(elem, circuit.vsource) or isinstance(elem, circuit.evsource) or \
		isinstance(elem, circuit.hvsource):
			#notice that hvsources aren't yet implemented now!
			i_eq = i_eq + 1
		elif isinstance(elem, circuit.capacitor):
			n1 = elem.n1
			n2 = elem.n2
			D[n1, n1] = D[n1, n1] + elem.C
			D[n1, n2] = D[n1, n2] - elem.C
			D[n2, n2] = D[n2, n2] + elem.C
			D[n2, n1] = D[n2, n1] - elem.C
		elif isinstance(elem, circuit.inductor):
			D[ nv + i_eq, nv + i_eq ] = -1 * elem.L
			i_eq = i_eq + 1
	return D

class dfbuffer:
	"""This is a LIFO buffer with a method to read it all without deleting the elements.
	Newer entries are added on top of the buffer.
	It checks the size of the added elements, to be sure they are of the same size.
	"""
	_the_real_buffer = None
	_length = 0
	_width  = 0
	
	def __init__(self, length, width):
		self._length = length
		self._width = width
		self._the_real_buffer = []
	
	def add(self, atuple):
		if not len(atuple) == self._width:
			printing.print_warning("Attempted to add a element of wrong size to LIFO buffer. BUG?")
			return False
		else:
			self._the_real_buffer.insert(0, atuple)
			if len(self._the_real_buffer) > self._length:
				self._the_real_buffer = self._the_real_buffer[:self._length]
			return True
	
	def get_df_vector(self):
		"""Returns a vector conforming to the specification of the df formulae. 
		That is [[time(n), x(n), dx(n)], [time(n-1), x(n-1), dx(n-1)], ...]
		"""
		return self._the_real_buffer
	
	def isready(self):
		"""This shouldn't be used to determine if the buffer has enough points to 
		use the df _if_ you use the step control.
		In that case, it holds even the points required for the FF.
		"""
		if len(self._the_real_buffer) == self._length:
			return True
		else:
			return False

def import_custom_df_module(method, print_out):
	"""Imports a module that implements differentiation formula through imp.load_module
	Parameters:
	method: a string, the name of the df method module
	print_out: print to stdout some verbose messages
	
	Returns:
	The df module or None if the module is not found.
	"""
	try:
		df = imp.load_module(imp.find_module(method.lower()))
		if print_out:
			print "Custom df module "+method.lower()+" loaded."
	except:
		printing.print_general_error("Unrecognized method: "+method.lower()+".")
		df = None
	
	return df
