# -*- coding: iso-8859-1 -*-
# dc_guess.py
# DC guess functions
# Copyright 2007 Giuseppe Venturini

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

"""This module provides the get_dc_guess() method.
"""

import sys
import numpy, numpy.linalg
import circuit, utilities

#try:
	#from scipy import factorial
#except ImportError:
	#from utilities import fact as factorial

def get_dc_guess(circ, verbose=3):
	"""This method tries to build a DC guess, according to what the
	elements suggest.
	A element can suggest its guess through the elem.dc_guess field.
	
	find_all: boolean, if set to True, the method will attempt to find all valid 
	guesses, it will not stop when it finds the first one.
	verbose: verbosity level (from 0 silent to 5 debug)
	
	Returns: a list of numpy's matrix
	"""
	
	#_debug = True

	
	if verbose: 
		sys.stdout.write("Calculating guess: ")
		sys.stdout.flush()
	
	# A DC guess has meaning only if the circuit has NL elements
	if not circ.is_nonlinear():
		if verbose:
			print "skipped. (linear circuit)"
		return []
	
	
	if verbose > 3:
		print ""

	nv = len(circ.nodes_dict)
	M = numpy.mat(numpy.zeros((1, nv)))
	T = numpy.mat(numpy.zeros((1, 1)))
	index = 0
	v_eq = 0 # number of current equations
	one_element_with_dc_guess_found = False

	for elem in circ.elements:
		# In the meanwhile, check how many current equations are 
		# required to solve the circuit
		if circuit.is_elem_voltage_defined(elem):
			v_eq = v_eq + 1
		# This is the main focus: build a system of equations (M*x = T)
		if hasattr(elem, "dc_guess") and elem.dc_guess is not None:
			if not one_element_with_dc_guess_found:
				one_element_with_dc_guess_found = True
			if elem.is_nonlinear:
				port_index = 0
				for (n1, n2) in elem.ports:
					if n1 == n2:
						continue
					if index:
						M = utilities.expand_matrix(M, add_a_row=True, add_a_col=False)
						T = utilities.expand_matrix(T, add_a_row=True, add_a_col=False)
					M[index, n1] = +1
					M[index, n2] = -1
					T[index] = elem.dc_guess[port_index]
					port_index = port_index + 1
					index = index + 1
			else:
				if elem.n1 == elem.n2:
					continue
				if index:
					M = utilities.expand_matrix(M, add_a_row=True, add_a_col=False)
					T = utilities.expand_matrix(T, add_a_row=True, add_a_col=False)
				M[index, elem.n1] = +1
				M[index, elem.n2] = -1
				T[index] = elem.dc_guess[0]
				index = index + 1
	
	if verbose == 5:
		print "DBG: get_dc_guess(): M and T, no reduction"
		print M
		print T
	M = utilities.remove_row_and_col(M, rrow=10*M.shape[0], rcol=0)
	
	if not one_element_with_dc_guess_found:
		if verbose == 5:
			print "DBG: get_dc_guess(): no element has a dc_guess"
		elif verbose <= 3:
			print "skipped."
		return []
	
	# We wish to find the linearly dependent lines of the M matrix.
	# The matrix is made by +1, -1, 0 elements. 
	# Hence, if two lines are linearly dependent, one of these equations
	# has to be satisfied: (L1, L2 are two lines)
	# L1 + L2 = 0 (vector)
	# L2 - L1 = 0 (vector)
	# This is tricky, because I wish to remove lines of the matrix while
	# browsing it.
	# We browse the matrix by line from bottom up and compare each line 
	# with the upper lines. If a linearly dep. line is found, we remove 
	# the current line.
	# Then break from the loop, get the next line (bottom up), which is
	# the same we were considering before; compare with the upper lines..
	# Not optimal, but it works.
	for i in range(M.shape[0]-1, -1, -1):
		#if not M[i, :].any(): #not needed -> does never happen
		#	M = utilities.remove_row(M, rrow=i)
		#	T = utilities.remove_row(T, rrow=i)
		#else:
		for j in range(i-1, -1, -1):
			#print i, j, M[i, :], M[j, :]
			dummy1 = M[i, :] - M[j, :]
			dummy2 = M[i, :] + M[j, :]
			if not dummy1.any() or not dummy2.any():
				#print "REM:", M[i, :]
				M = utilities.remove_row(M, rrow=i)
				T = utilities.remove_row(T, rrow=i)
				break
	if verbose == 5:
		print "DBG: get_dc_guess(): M and T, after removing LD lines"
		print M
		print T
			
	# Remove empty columns:
	# If a column is empty, we have no guess regarding the corresponding
	# node. It makes the matrix singular. -> Remove the col & remember
	# that we are _not_ calculating a guess for it.
	removed_index = []
	for i in range(M.shape[1]-1, -1, -1):
		if not M[:, i].any():
			M = utilities.remove_row_and_col(M, rrow=M.shape[0], rcol=i)
			removed_index.append(i)
	
	if verbose > 3:
		print "DBG: get_dc_guess(): M and T, after removing empty columns."
		print M
		print "T\n", T

	# Now, we have a set of equations to be solved.
	# There are three cases:
	# 1. The M matrix has more rows than columns
	#    If this happens, it means that there is no solution that satisfies
	#    all the guesses correctly.
	#    We have to choose a subset of the rows them to get a square M. 
	#    Notice that:
	#     -> Different choices _may_ give different solutions.
	#     -> Some choices may give a singular M matrix.
	# 2. The matrix has more columns than rows.
	#    This means that we haven't got enough information to solve for all nodes.
	#    We have to choose a subset of the columns (and of the variables).
	#    The others will be left to 0.
	# 3. The matrix is square.
	#    I'm not sure about this: it seems that if the circuit is not pathological,
	#    we are likely to find a solution (the matrix has det != 0).
	
	Rp_list = [] # to hold the results
	if M.shape[0] != M.shape[1]:
		Rp_list = [numpy.mat(numpy.linalg.pinv(M))*T]
	else: # case M.shape[0] == M.shape[1], use normal
		if numpy.linalg.det(M) is not 0:
			try:
				Rp_list = Rp_list + [numpy.linalg.inv(M) * T]
			except numpy.linalg.linalg.LinAlgError:
				eig = numpy.linalg.eig(M)[0]
				cond = abs(eig).max()/abs(eig).min()
				if verbose and verbose < 4:
					sys.stdout.write("cond=" +str(cond)+". No guess. ")

	# Now we want to:
	# 1. Add voltages for the nodes for which we have no clue to guess.
	# 2. Append to each vector of guesses the values for currents in 
	#    voltage defined elem.
	# Both them are set to 0
	Rp_list_dummy = [] # We want to modify and readd each item
	for Rp in Rp_list:
		for index in removed_index:
			Rp = numpy.concatenate(( \
			numpy.concatenate((Rp[:index, 0], \
			numpy.mat(numpy.zeros((1, 1)))), axis=0), \
			Rp[index:, 0]), axis=0)
		# add the 0s for the currents due to the voltage defined 
		# elements (we have no guess for those...)
		if v_eq > 0:
			Rp = numpy.concatenate((Rp, numpy.mat(numpy.zeros((v_eq, 1)))), axis=0)
		Rp_list_dummy = Rp_list_dummy + [Rp]
	Rp_list = Rp_list_dummy
	del Rp_list_dummy
	
	if verbose == 5:
		print circ.nodes_dict
		print Rp_list
	
	if verbose and verbose < 4:
		print "done."
	if verbose > 3:
		print "Guesses:", len(Rp_list)
		for Rp in Rp_list:
			print Rp
	
	return Rp_list
