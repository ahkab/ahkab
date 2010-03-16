import sympy
import circuit, ekv, mosq, printing

def dc_solve(circ, ac=False, tf_source=None, options={'r0s':True}, verbose=6):
	#options = setup_options()	
	if verbose > 1 and not ac:
		print "Starting symbolic DC..."
	elif verbose > 1 and ac:
		print "Starting symbolic AC..."
		
	if verbose > 2:
		 print "Building symbolic MNA, N and x..."
	mna, N = generate_mna_and_N(circ, options)
	x = generate_variable_names(circ, N.shape[0] - 1)
	mna = mna[1:, 1:]
	N = N[1:, :]
	
	#if verbose > 2:
	#	 print "Building equation..."
	#x = to_real_list(x)
	#eq = apply_options(eq, options)

	if verbose > 2:
		 print "Solving..."
	if verbose > 5:
		eq = to_real_list(mna * x + N)
		printing.print_symbolic_equations(eq)
		print "To be solved for:"
		print to_real_list(x)
		#print mna.inv()
		print "Matrix is singular: ", (mna.det() == 0)
	sol = -1.0*mna.inv()*N
	sol = sol_to_dict(sol, x)
	#sol = sympy.solve(eq, x)
	
	if sol == {}:
		print "HEY, NO VARIABLES?"
	else:
		if verbose > 2:
		 	print "Success!"
		elif verbose > 1:
			print "Results:"
		printing.print_symbolic_results(sol)
	
	if tf_source is not None:
		if verbose > 2: print "Calculating symbolic transfer functions...",
		src = sympy.Symbol(tf_source, real=True)
		tfs = calculate_gains(sol, sympy.Symbol("v2", real=True))
		if verbose > 2: print "done!"
		elif verbose > 1:
			print "Symbolic transfer functions:"
		printing.print_symbolic_results(tfs)
	

def calculate_gains(sol, xin, optimize=True):
	gains = {}
	for key, value in sol.iteritems():
		gain = sympy.simplify(value.diff(xin)) if optimize else value.diff(xin)
		gains.update({"d/d" + str(xin)+" "+str(key):gain})
	return gains


def sol_to_dict(sol, x, optimize=False):
	ret = {}
	for index in range(x.shape[0]):
		sol_current = sympy.simplify(sol[index]) if optimize else sol[index]
		ret.update({str(x[index]):sol_current})
	return ret

def apply_options(eq_list, options):
	running_eq_list = eq_list	
	for key in options.keys():
		subs_eq_list = []		
		for eq in running_eq_list:
			subs_eq_list.append(eq.subs(key, options[key]))
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
	options = {}	
	return options

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

def generate_mna_and_N(circ, options, ac=False):
	"""Generates a symbolic Modified Nodal Analysis matrix and the N vector.
	"""
	#print options
	n_of_nodes = len(circ.nodes_dict)
	mna = sympy.matrices.zeros(n_of_nodes)
	N = sympy.matrices.zeros((n_of_nodes, 1))
	omega = sympy.Symbol("w", real=True)
	#process_elements() 	
	for elem in circ.elements:
		if elem.is_nonlinear and not (isinstance(elem, mosq.mosq) or isinstance(elem, ekv.ekv_device)): 
			print "Skipped elem "+elem.letter_id+elem.descr	
			continue
		if isinstance(elem, circuit.resistor):
			R = sympy.Symbol(elem.letter_id+elem.descr)
			mna[elem.n1, elem.n1] = mna[elem.n1, elem.n1] + 1/R
			mna[elem.n1, elem.n2] = mna[elem.n1, elem.n2] - 1/R
			mna[elem.n2, elem.n1] = mna[elem.n2, elem.n1] - 1/R
			mna[elem.n2, elem.n2] = mna[elem.n2, elem.n2] + 1/R
		elif isinstance(elem, circuit.capacitor) and ac:
			capa = sympy.Symbol(elem.letter_id+elem.descr, real=True)
			mna[elem.n1, elem.n1] = mna[elem.n1, elem.n1] + sympy.I*omega*capa
			mna[elem.n1, elem.n2] = mna[elem.n1, elem.n2] - sympy.I*omega*capa
			mna[elem.n2, elem.n2] = mna[elem.n2, elem.n2] + sympy.I*omega*capa
			mna[elem.n2, elem.n1] = mna[elem.n2, elem.n1] - sympy.I*omega*capa
		elif isinstance(elem, circuit.inductor):
			pass
		elif isinstance(elem, circuit.gisource):
			alpha = sympy.Symbol(elem.letter_id+elem.descr, real=True)
			mna[elem.n1, elem.sn1] = mna[elem.n1, elem.sn1] + alpha
			mna[elem.n1, elem.sn2] = mna[elem.n1, elem.sn2] - alpha
			mna[elem.n2, elem.sn1] = mna[elem.n2, elem.sn1] - alpha
			mna[elem.n2, elem.sn2] = mna[elem.n2, elem.sn2] + alpha
		elif isinstance(elem, circuit.isource):
			IDC = sympy.Symbol(elem.letter_id+elem.descr, real=True)
			N[elem.n1, 0] = N[elem.n1, 0] + IDC
			N[elem.n2, 0] = N[elem.n2, 0] - IDC
		elif isinstance(elem, mosq.mosq) or isinstance(elem, ekv.ekv_device):
			gm = sympy.Symbol('gm_'+elem.letter_id+elem.descr, real=True)
			mna[elem.n1, elem.ng] = mna[elem.n1, elem.ng] + gm
			mna[elem.n1, elem.n2] = mna[elem.n1, elem.n2] - gm
			mna[elem.n2, elem.ng] = mna[elem.n2, elem.ng] - gm
			mna[elem.n2, elem.n2] = mna[elem.n2, elem.n2] + gm
			if options['r0s']:
				r0 = sympy.Symbol('r0_'+elem.letter_id+elem.descr, real=True)
				mna[elem.n1, elem.n1] = mna[elem.n1, elem.n1] + 1/r0
				mna[elem.n1, elem.n2] = mna[elem.n1, elem.n2] - 1/r0
				mna[elem.n2, elem.n1] = mna[elem.n2, elem.n1] - 1/r0
				mna[elem.n2, elem.n2] = mna[elem.n2, elem.n2] + 1/r0
		elif circuit.is_elem_voltage_defined(elem):
			pass
			#we'll add its lines afterwards
		else:
			print "dcsymbolic.py: BUG - Unknown linear element. Ref. #28935"

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
			if isinstance(elem, circuit.vsource):
				N[index, 0] = -1.0 * sympy.Symbol(elem.letter_id + elem.descr, real=True)
			elif isinstance(elem, circuit.evsource):
				alpha = sympy.Symbol(elem.letter_id + elem.descr, real=True)
				mna[index, elem.sn1] = -1.0 * alpha
				mna[index, elem.sn2] = +1.0 * alpha
			elif isinstance(elem, circuit.inductor) and ac:
				mna[index, index] = -1*sympy.I*omega* sympy.Symbol(elem.letter_id + elem.descr, real=True)
				# already so: commented out				
				# N[index,0] = 0
			elif isinstance(elem, circuit.hvsource):
				print "dcsymbolic.py: BUG - hvsources are not implemented yet."
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

#def process_elements(circ):
#	new_elem_list = []
#	for elem in circ.elements:		
#		if isinstance(elem, mosq.mosq):
#			circuit.resistor(elem.nd, elem.ns)			
#			else:

#def build_mos_function(vg, vs, )
