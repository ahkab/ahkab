# -*- coding: iso-8859-1 -*-
# symbolic.py
# Symbolic simulation module
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

"""
This module offers the functions needed to perform a symbolic simulation,
AC or DC.

The principal is solve() - which carries out the symbolic solution
"""

import sympy
import circuit, devices, ekv, mosq, printing, options

def solve(circ, ac=False, tf_source=None, opts={'r0s':True}, verbose=6):
	#opts = setup_options()	
	if verbose > 1 and not ac:
		print "Starting symbolic DC..."
	elif verbose > 1 and ac:
		print "Starting symbolic AC..."
		
	if verbose > 2:
		 print "Building symbolic MNA, N and x..."
	mna, N = generate_mna_and_N(circ, opts, ac)
	x = generate_variable_names(circ, N.shape[0] - 1)
	mna = mna[1:, 1:]
	N = N[1:, :]
	
	if verbose > 5:
		print "MNA matrix (reduced):"
		print mna
	
	if verbose > 2:
		 print "Building equation..."
	eq = to_real_list(mna * x + N)
	#eq = apply_options_and_subs(eq, opts)

	x = to_real_list(x)
	if verbose > 4:
		printing.print_symbolic_equations(eq)
		print "To be solved for:"
		print x
		#print "Matrix is singular: ", (mna.det() == 0)
	#sol = -1.0*mna.inv()*N #too heavy
	if verbose > 2:
		print "Performing auxiliary simplification..."
	eq, x, sol_h = help_the_solver(eq, x)
		
	if len(eq):
		if verbose > 3:
			print "Symplified sytem:"
			printing.print_symbolic_equations(eq)
			print "To be solved for:"
			print x
		if verbose > 2:
			print "Solving..."

		if options.symb_internal_solver:
			sol = local_solve(eq, x)
		else:
			sol = sympy.solve(eq, x)
		sol.update(sol_h)
	else:
		if verbose > 3:
			print "Auxiliary simplification solved the problem."
		sol = sol_h
	
	if verbose > 2:
	 	print "Success!"
	#sol = sol_to_dict(sol, x)
	
	if sol == {}:
		print "HEY, NO VARIABLES?"
	else:
		if verbose > 1:
			print "Results:"
		printing.print_symbolic_results(sol)
	
	if tf_source is not None:
		src = sympy.Symbol(tf_source, real=True)
		if verbose > 2: print "Calculating small-signal symbolic transfer functions ("+str(src)+")...",
		tfs = calculate_gains(sol, src)
		if verbose > 2: print "done!"
		elif verbose > 1:
			print "Small-signal symbolic transfer functions:"
		printing.print_symbolic_transfer_functions(tfs)
	

def calculate_gains(sol, xin, optimize=True):
	gains = {}
	for key, value in sol.iteritems():
		tf = {}
		gain = sympy.together(value.diff(xin)) if optimize else value.diff(xin)
		(ps, zs) = get_roots(gain)
		tf.update({'gain':gain})
		tf.update({'gain0':gain.subs(sympy.Symbol('s', real=False), 0)})
		tf.update({'poles':ps})
		tf.update({'zeros':zs})
		gains.update({"%s/%s" % (str(key), str(xin)):tf})
	return gains


def sol_to_dict(sol, x, optimize=True):
	ret = {}
	for index in range(x.shape[0]):
		sol_current = sympy.together(sol[index]) if optimize else sol[index]
		ret.update({str(x[index]):sol_current})
	return ret

def apply_options_and_subs(eq_list, opts):
	import subs
	opts = {}
	for key, value in subs.subs.iteritems():
		opts.update({sympy.Symbol(key, real=True):sympy.Symbol(value, real=True)}) 
	subs_eq_list = []		
	running_eq_list = eq_list	
	for key in opts.keys():
		subs_eq_list = []		
		for eq in running_eq_list:
			subs_eq_list.append(eq.subs(key, opts[key]))
		running_eq_list = subs_eq_list
	print subs_eq_list
	return subs_eq_list

def setup_options():
	"""options = {}	
	r1 = sympy.Symbol('r1')
	r2 = sympy.Symbol('r2')
	r3 = sympy.Symbol('r3')
	r4 = sympy.Symbol('r4')
	r5 = sympy.Symbol('r5')
	r6 = sympy.Symbol('r6')
	r7 = sympy.Symbol('r7')
	r8 = sympy.Symbol('r8')
	r9 = sympy.Symbol('r9')
	r10 = sympy.Symbol('r10')
	options.update({r2:2*r1})
	options.update({r3:r1})
	options.update({r4:2*r1})
	options.update({r5:r1})
	options.update({r6:2*r1})
	options.update({r7:r1})
	options.update({r8:2*r1})
	options.update({r9:r1})
	options.update({r10:r1})"""
	opts = {}	
	return opts

def generate_variable_names(circ, mna_size):
	x = sympy.matrices.zeros((mna_size, 1))

	nv_1 = len(circ.nodes_dict) - 1 # numero di soluzioni di tensione (al netto del ref)
	
	# descrizioni dei componenti non definibili in tensione
	idescr = [ (elem.letter_id.upper() + elem.descr) \
		for elem in circ.elements if circuit.is_elem_voltage_defined(elem) ]

	for i in range(mna_size):

		if i < nv_1:
			x[i, 0] = sympy.Symbol("V" + str(circ.nodes_dict[i + 1]))
		else:
			x[i, 0] = sympy.Symbol("I["+idescr[i - nv_1]+"]")
	return x

def to_real_list(M):
	"""
	M.tolist() returns a list of lists, even when the symb matrix is really just a vector.
	we want a list of symbols! This fixes that.

	mylist[k] = mymat.tolist[k][0]

	M: a sympy matrix with only one column

	Returns: a list.
	"""
	fakelist = M.tolist();
	reallist = []
	for elem in fakelist:
		reallist.append(elem[0])
	return reallist

def generate_mna_and_N(circ, opts, ac=False):
	"""Generates a symbolic Modified Nodal Analysis matrix and the N vector.
	"""
	#print options
	n_of_nodes = len(circ.nodes_dict)
	mna = sympy.matrices.zeros(n_of_nodes)
	N = sympy.matrices.zeros((n_of_nodes, 1))
	s = sympy.Symbol("s", real=False)
	#process_elements() 	
	for elem in circ.elements:
		#if elem.is_nonlinear and not (isinstance(elem, mosq.mosq) or isinstance(elem, ekv.ekv_device)): 
		#	print "Skipped elem "+elem.letter_id+elem.descr + ": not implemented."	
		#	continue
		if isinstance(elem, devices.resistor):
			R = sympy.Symbol(elem.letter_id.upper()+elem.descr)
			mna[elem.n1, elem.n1] = mna[elem.n1, elem.n1] + 1/R
			mna[elem.n1, elem.n2] = mna[elem.n1, elem.n2] - 1/R
			mna[elem.n2, elem.n1] = mna[elem.n2, elem.n1] - 1/R
			mna[elem.n2, elem.n2] = mna[elem.n2, elem.n2] + 1/R
		elif isinstance(elem, devices.capacitor):
			if ac:
				capa = sympy.Symbol(elem.letter_id.upper()+elem.descr, real=True)
				mna[elem.n1, elem.n1] = mna[elem.n1, elem.n1] + s*capa
				mna[elem.n1, elem.n2] = mna[elem.n1, elem.n2] - s*capa
				mna[elem.n2, elem.n2] = mna[elem.n2, elem.n2] + s*capa
				mna[elem.n2, elem.n1] = mna[elem.n2, elem.n1] - s*capa
			else:
				pass
		elif isinstance(elem, devices.inductor):
			pass
		elif isinstance(elem, devices.gisource):
			alpha = sympy.Symbol(elem.letter_id+elem.descr, real=True)
			mna[elem.n1, elem.sn1] = mna[elem.n1, elem.sn1] + alpha
			mna[elem.n1, elem.sn2] = mna[elem.n1, elem.sn2] - alpha
			mna[elem.n2, elem.sn1] = mna[elem.n2, elem.sn1] - alpha
			mna[elem.n2, elem.sn2] = mna[elem.n2, elem.sn2] + alpha
		elif isinstance(elem, devices.isource):
			IDC = sympy.Symbol(elem.letter_id+elem.descr, real=True)
			N[elem.n1, 0] = N[elem.n1, 0] + IDC
			N[elem.n2, 0] = N[elem.n2, 0] - IDC
		elif isinstance(elem, mosq.mosq) or isinstance(elem, ekv.ekv_device):
			gm = sympy.Symbol('gm_'+elem.letter_id+elem.descr, real=True)
			mna[elem.n1, elem.ng] = mna[elem.n1, elem.ng] + gm
			mna[elem.n1, elem.n2] = mna[elem.n1, elem.n2] - gm
			mna[elem.n2, elem.ng] = mna[elem.n2, elem.ng] - gm
			mna[elem.n2, elem.n2] = mna[elem.n2, elem.n2] + gm
			if opts['r0s']:
				r0 = sympy.Symbol('r0_'+elem.letter_id+elem.descr, real=True)
				mna[elem.n1, elem.n1] = mna[elem.n1, elem.n1] + 1/r0
				mna[elem.n1, elem.n2] = mna[elem.n1, elem.n2] - 1/r0
				mna[elem.n2, elem.n1] = mna[elem.n2, elem.n1] - 1/r0
				mna[elem.n2, elem.n2] = mna[elem.n2, elem.n2] + 1/r0
		elif isinstance(elem, devices.diode):
			gd = sympy.Symbol("g"+elem.letter_id+elem.descr)
			mna[elem.n1, elem.n1] = mna[elem.n1, elem.n1] + gd
			mna[elem.n1, elem.n2] = mna[elem.n1, elem.n2] - gd
			mna[elem.n2, elem.n1] = mna[elem.n2, elem.n1] - gd
			mna[elem.n2, elem.n2] = mna[elem.n2, elem.n2] + gd
		elif circuit.is_elem_voltage_defined(elem):
			pass
			#we'll add its lines afterwards
		else:
			print "Skipped elem "+elem.letter_id+elem.descr + ": not implemented."	

	for elem in circ.elements:
		if circuit.is_elem_voltage_defined(elem):
			index = mna.shape[0] #get_matrix_size(mna)[0]
			mna = expand_matrix(mna, add_a_row=True, add_a_col=True)
			N = expand_matrix(N, add_a_row=True, add_a_col=False)
			# KCL
			mna[elem.n1, index] = +1.0
			mna[elem.n2, index] = -1.0
			# KVL
			mna[index, elem.n1] = +1.0
			mna[index, elem.n2] = -1.0
			if isinstance(elem, devices.vsource):
				N[index, 0] = -1.0 * sympy.Symbol(elem.letter_id + elem.descr, real=True)
			elif isinstance(elem, devices.evsource):
				alpha = sympy.Symbol(elem.letter_id + elem.descr, real=True)
				mna[index, elem.sn1] = -1.0 * alpha
				mna[index, elem.sn2] = +1.0 * alpha
			elif isinstance(elem, devices.inductor):
				if ac:
					mna[index, index] = -1*s*sympy.Symbol(elem.letter_id.upper() + elem.descr, real=True)
				else: 
					pass
					# already so: commented out				
					# N[index,0] = 0
			elif isinstance(elem, devices.hvsource):
				print "symbolic.py: BUG - hvsources are not implemented yet."
				sys.exit(33)
	#all done
	return (mna, N)

def expand_matrix(mat, add_a_row=False, add_a_col=False):
	if add_a_row:
		row = sympy.zeros((1, mat.shape[1]))
		mat = mat.row_insert(mat.shape[0], row)
	if add_a_col:
		col = sympy.zeros((mat.shape[0], 1))
		mat = mat.col_insert(mat.shape[1], col)
	return mat

def get_roots(expr):
	num, den = sympy.fraction(expr)
	s = sympy.Symbol('s', real=False)
	return sympy.solve(den, s), sympy.solve(num, s)

############## THESE  WILL BE REMOVED - AS SOON AS SOME SYMPY BUGS ARE FIXED ###########
def help_the_solver(eqs, xs, debug=True):
	iter_flag = True
	sol = {}
	while iter_flag:
		iter_flag, eqs, subs = help_the_solver_iter(eqs, xs)
		if iter_flag:
			xs.remove(subs.keys()[0])
			sol.update(subs)
		if debug:
			for key, value in subs.iteritems():
				print key, "=", value
	return eqs, xs, sol

def help_the_solver_iter(eqs, xs):
	success = False
	for eq in eqs:
		success, subs = help_the_solver_1eq(eq, xs)
		if success:
			break
	if success:
		new_eqs = []
		eqs.remove(eq)
		for eq in eqs:
			new_eqs.append(eq.subs(subs))
	else:
		new_eqs = eqs
		subs = {}
	return success, new_eqs, subs
def help_the_solver_1eq(eq, xs, debug=True):
	one_x = None
	for x in xs:
		if eq.has(x) and one_x is None:
			one_x = x
		elif eq.has(x) and one_x is not None:
			one_x = None
			break
	if one_x is not None:
		sol = {one_x:sympy.solve(eq, one_x)[0]}
	else:
		sol = {}	
	return not one_x is None, sol

def local_solve(eqs, xs):
	sol = {}
	while len(eqs):
		eqs, single_sol = local_solve_iter(eqs, xs)
		new_sol = {}
		for key, value in sol.iteritems():
			new_sol.update({key:value.subs(single_sol)})
		new_sol.update(single_sol)
		sol = new_sol
		new_eqs = []
		for eq in eqs:
			new_eqs.append(eq.subs(single_sol)) 
		eqs = new_eqs
	return sol

def local_solve_iter(eqs, xs):
	for eq in eqs:	
		for x in xs:
			if eq.has(x):
				print "Solving for", x
				single_sol = {x:sympy.solve(eq, x)[0]}
				eqs.remove(eq)
				print single_sol
				return eqs, single_sol
	return eqs, {}
	

#def process_elements(circ):
#	new_elem_list = []
#	for elem in circ.elements:		
#		if isinstance(elem, mosq.mosq):
#			devices.resistor(elem.nd, elem.ns)			
#			else:

#def build_mos_function(vg, vs, )
