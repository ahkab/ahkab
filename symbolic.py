# -*- coding: iso-8859-1 -*-
# symbolic.py
# Symbolic simulation module
# Copyright 2010-2013 Giuseppe Venturini

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

The principal method is solve(), which carries out the symbolic circuit solution.
"""

import sympy
from sympy.matrices import zeros as smzeros

import circuit
import devices 
import ekv 
import mosq
import printing
import options

def solve(circ, tf_source=None, subs=None, opts=None, verbose=3):
	"""Attempt a symbolic solution of the circuit.
	circ: the circuit instance to be simulated.
	tf_source: the name (string) of the source to be used as input for the transfer
		   function. If None, no transfer function is evaluated.
	subs: a dictionary of sympy Symbols to be substituted. It makes solving the circuit 
	      easier. Eg. {R1:R2} - replace R1 with R2. It can be generated with 
	      parse_substitutions()
	opts: dict of 'option':boolean to be taken into account in simulation.
	      currently 'r0s' and 'ac' are the only options considered.
	verbose: verbosity level 0 (silent) to 6 (painful).
	
	Returns: a dictionary with the solutions.
	"""
	if opts is None:
		# load the defaults
		opts = {'r0s':True, 'ac':False}
	if not 'r0s' in opts.keys():
		opts.update({'r0s':True})
	if not 'ac' in opts.keys():
		opts.update({'ac':False})
	if subs is None:
		subs = {} # no subs by default

	if not opts['ac']:
		printing.print_info_line(("Starting symbolic DC...", 1), verbose)
	else:
		printing.print_info_line(("Starting symbolic AC...", 1), verbose)		
		
	printing.print_info_line(("Building symbolic MNA, N and x...", 2), verbose, print_nl=False)
	mna, N, subs_g = generate_mna_and_N(circ, opts, opts['ac'])
	x = get_variables(circ)
	mna = mna[1:, 1:]
	N = N[1:, :]
	printing.print_info_line((" done.", 2), verbose)	

	printing.print_info_line(("Performing variable substitutions...", 5), verbose)
	mna, N = apply_substitutions(mna, N, subs)

	printing.print_info_line(("MNA matrix (reduced):", 5), verbose)	
	if verbose > 5:	print sympy.sstr(mna)

	printing.print_info_line(("Building equations...", 2), verbose)	
	eq = to_real_list(mna * x + N)

	x = to_real_list(x)
	if verbose > 4:
		printing.print_symbolic_equations(eq)
		print "To be solved for:"
		print x
		#print "Matrix is singular: ", (mna.det() == 0)
	#print -1.0*mna.inv()*N #too heavy
	#print sympy.solve_linear_system(mna, x)
	printing.print_info_line(("Performing auxiliary simplification...", 2), verbose)	
	eq, x, sol_h = help_the_solver(eq, x)
		
	if len(eq):
		if verbose > 3:
			print "Symplified sytem:"
			printing.print_symbolic_equations(eq)
			print "To be solved for:"
			print x
			printing.print_info_line(("Solving...", 2), verbose)	

		if options.symb_internal_solver:
			sol = local_solve(eq, x)
		else:
			sol = sympy.solve(eq, x, manual=options.symb_sympy_manual_solver, simplify=False)
		if sol is not None:		
			sol.update(sol_h)
		else:
			sol = sol_h
	else:
		printing.print_info_line(("Auxiliary simplification solved the problem.", 3), verbose)	
		sol = sol_h

	for ks in sol.keys():
		sol.update({ks:sol[ks].subs(subs_g)})

	#sol = sol_to_dict(sol, x)

	if sol == {}:
		printing.print_warning("No solutions. Check the netlist.")
	else:
		printing.print_info_line(("Success!", 2), verbose)	
		if verbose > 1:
			print "Results:"
		printing.print_symbolic_results(sol)

	if tf_source is not None:
		src = sympy.Symbol(tf_source.lower(), real=True)
		printing.print_info_line(("Calculating small-signal symbolic transfer functions (%s))..."%(str(src),), 2), verbose, print_nl=False)
		tfs = calculate_gains(sol, src)
		printing.print_info_line(("done.", 2), verbose)	
		printing.print_info_line(("Small-signal symbolic transfer functions:", 1), verbose)	
		printing.print_symbolic_transfer_functions(tfs)
	else:
		tfs = None
	
	return sol, tfs

def calculate_gains(sol, xin, optimize=True):
	gains = {}
	for key, value in sol.iteritems():
		tf = {}
		gain = sympy.together(value.diff(xin)) if optimize else value.diff(xin)
		(ps, zs) = get_roots(gain)
		tf.update({'gain':gain})
		tf.update({'gain0':gain.subs(sympy.Symbol('s', complex=True), 0)})
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

def apply_substitutions(mna, N, opts):
	mna = mna.subs(opts)
	N = N.subs(opts)
	return (mna, N)

def get_variables(circ):
	"""Returns a list with the circuit variables to be solved for.
	"""
	nv_1 = len(circ.nodes_dict) - 1 # numero di soluzioni di tensione (al netto del ref)
	
	# descrizioni dei componenti non definibili in tensione
	idescr = [ (elem.letter_id.upper() + elem.descr) \
		for elem in circ.elements if circuit.is_elem_voltage_defined(elem) ]

	mna_size = nv_1 + len(idescr)
	x = smzeros((mna_size, 1))

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
	"""Generates a symbolic Modified Nodal Analysis matrix and N vector.
	"""
	#print options
	n_of_nodes = len(circ.nodes_dict)
	mna = smzeros(n_of_nodes)
	N = smzeros((n_of_nodes, 1))
	s = sympy.Symbol("s", complex=True)
	subs_g = {}
	#process_elements() 	
	for elem in circ.elements:
		#if elem.is_nonlinear and not (isinstance(elem, mosq.mosq_device) or isinstance(elem, ekv.ekv_device)): 
		#	print "Skipped elem "+elem.letter_id+elem.descr + ": not implemented."	
		#	continue
		if isinstance(elem, devices.resistor):
			# we use conductances instead of 1/R because there is a significant
			# overhead handling many 1/R terms in sympy.
			R = sympy.Symbol(elem.letter_id.upper()+elem.descr, real=True)
			#mna[elem.n1, elem.n1] = mna[elem.n1, elem.n1] + 1/R
			#mna[elem.n1, elem.n2] = mna[elem.n1, elem.n2] - 1/R
			#mna[elem.n2, elem.n1] = mna[elem.n2, elem.n1] - 1/R
			#mna[elem.n2, elem.n2] = mna[elem.n2, elem.n2] + 1/R
			G = sympy.Symbol("G"+elem.descr, real=True)
			mna[elem.n1, elem.n1] = mna[elem.n1, elem.n1] + G
			mna[elem.n1, elem.n2] = mna[elem.n1, elem.n2] - G
			mna[elem.n2, elem.n1] = mna[elem.n2, elem.n1] - G
			mna[elem.n2, elem.n2] = mna[elem.n2, elem.n2] + G
			# but we keep track of which is which and substitute back after solving.
			subs_g.update({G:1/R})
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
		elif isinstance(elem, mosq.mosq_device) or isinstance(elem, ekv.ekv_device):
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
		elif isinstance(elem, devices.inductor_coupling):
			pass
			# this is taken care of within the inductors
		elif circuit.is_elem_voltage_defined(elem):
			pass
			#we'll add its lines afterwards
		else:
			printing.print_warning("Skipped elem %s: not implemented." % (elem.letter_id+elem.descr,))

	pre_vde = mna.shape[0]
	for elem in circ.elements:
		if circuit.is_elem_voltage_defined(elem):
			index = mna.shape[0] #get_matrix_size(mna)[0]
			mna = expand_matrix(mna, add_a_row=True, add_a_col=True)
			N = expand_matrix(N, add_a_row=True, add_a_col=False)
			# KCL
			mna[elem.n1, index] = +1
			mna[elem.n2, index] = -1
			# KVL
			mna[index, elem.n1] = +1
			mna[index, elem.n2] = -1
			if isinstance(elem, devices.vsource):
				N[index, 0] = -sympy.Symbol(elem.letter_id + elem.descr, real=True)
			elif isinstance(elem, devices.evsource):
				alpha = sympy.Symbol(elem.letter_id.upper() + elem.descr, real=True)
				mna[index, elem.sn1] = -alpha
				mna[index, elem.sn2] = +alpha
			elif isinstance(elem, devices.inductor):
				if ac:
					mna[index, index] = -1*s*sympy.Symbol(elem.letter_id.upper() + elem.descr, real=True)
				else: 
					pass
					# already so: commented out				
					# N[index,0] = 0
			elif isinstance(elem, devices.hvsource):
				printing.print_warning("symbolic.py: BUG - hvsources are not implemented yet.")
				sys.exit(33)
	
	for elem in circ.elements:
		if circuit.is_elem_voltage_defined(elem):
			if isinstance(elem, devices.inductor):
				if ac:
					# find its index to know which column corresponds to its current
					this_index = circ.find_vde_index("L"+elem.descr, verbose=0)
					for cd in elem.coupling_devices:
						# get id+descr of the other inductor (eg. "L32")
						other_id_wdescr = cd.get_other_inductor("L"+elem.descr)
						# find its index to know which column corresponds to its current
						other_index = circ.find_vde_index(other_id_wdescr, verbose=0)
						# add the term.
						#print "other_index: "+str(other_index)
						#print "this_index: "+str(this_index)
						mna[pre_vde+this_index,pre_vde+other_index] += -1*s*sympy.Symbol("M" + cd.descr, real=True)
						#print mna
				else: 
					pass
					# already so: commented out				
					# N[index,0] = 0
			
	
	#all done
	return (mna, N, subs_g)

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
	s = sympy.Symbol('s', complex=True)
	return sympy.solve(den, s), sympy.solve(num, s)

def parse_substitutions(slist):
	"""Generates a substitution dictionary to be passed to solve()
	
	slist is a list of strings like 'R2=R1', instructing the simulator
	to use the value of R1 (R1) instead of R2.

	returns: the dictionary of symbols to be passed to solve()
	"""
	subs = {}
	for l in slist:
		v1, v2 = l.split("=")
		letter_id1 = v1[0].upper() if v1[0].upper() != 'R' else 'G'
		letter_id2 = v2[0].upper() if v2[0].upper() != 'R' else 'G'
		subs.update({
			sympy.Symbol(letter_id1+v1[1:].lower(), real=True):
			sympy.Symbol(letter_id2+v2[1:].lower(), real=True)
			})
	return subs

############## THESE  WILL BE REMOVED - AS SOON AS SOME SYMPY BUGS ARE FIXED ###########
def help_the_solver(eqs, xs, debug=False):
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
