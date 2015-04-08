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

"""This module provides the :func:`get_dc_guess` method, used to
compute a starting point to initialize a Newton-Rhapson solver.

Module reference
################

"""

from __future__ import (unicode_literals, absolute_import,
                        division, print_function)

import sys

import numpy as np
import numpy.linalg

from . import circuit
from . import utilities


def get_dc_guess(circ, verbose=3):
    """Build a DC guess from circuit inspection.

    Notice that OP analysis will call this method on the users' behalf if not
    instructed not to do so.

    A element can suggest its guess through the ``elem.dc_guess`` field.
    If the field is not set, or not available, no information on the most
    likely biasing voltage is assumed.

    **Parameters:**

    circ : Circuit instance
        The circuit instance the guess is being computed for.

    verbose : int, optional
        The verbosity level (from 0 silent to 6 debug). Defaults to 3, medium
        verbosity.

    **Returns:**

    dcg : ndarray or None
        The DC guess, in numpy array form, or ``None``, if it was not possible
        to compute a meaningful guess.
    """
    if verbose:
        sys.stdout.write("Calculating guess: ")
        sys.stdout.flush()

    # A DC guess has meaning only if the circuit has NL elements
    if not circ.is_nonlinear():
        if verbose:
            print("skipped. (linear circuit)")
        return None

    if verbose > 3:
        print("")

    nv = circ.get_nodes_number()
    M = np.zeros((1, nv))
    T = np.zeros((1, 1))
    index = 0
    v_eq = 0  # number of current equations
    one_element_with_dc_guess_found = False

    for elem in circ:
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
                        M = utilities.expand_matrix(M, add_a_row=True,
                                                    add_a_col=False)
                        T = utilities.expand_matrix(T, add_a_row=True,
                                                    add_a_col=False)
                    M[index, n1] = +1
                    M[index, n2] = -1
                    T[index] = elem.dc_guess[port_index]
                    port_index = port_index + 1
                    index = index + 1
            else:
                if elem.n1 == elem.n2:
                    continue
                if index:
                    M = utilities.expand_matrix(M, add_a_row=True,
                                                add_a_col=False)
                    T = utilities.expand_matrix(T, add_a_row=True,
                                                add_a_col=False)
                M[index, elem.n1] = +1
                M[index, elem.n2] = -1
                T[index] = elem.dc_guess[0]
                index = index + 1

    if verbose == 5:
        print("DBG: get_dc_guess(): M and T, no reduction")
        print(M)
        print(T)
    M = utilities.remove_row_and_col(M, rrow=10 * M.shape[0], rcol=0)

    if not one_element_with_dc_guess_found:
        if verbose == 5:
            print("DBG: get_dc_guess(): no element has a dc_guess")
        elif verbose <= 3:
            print("skipped.")
        return None

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
    for i in range(M.shape[0] - 1, -1, -1):
        for j in range(i - 1, -1, -1):
            # print i, j, M[i, :], M[j, :]
            dummy1 = M[i, :] - M[j, :]
            dummy2 = M[i, :] + M[j, :]
            if not dummy1.any() or not dummy2.any():
                # print "REM:", M[i, :]
                M = utilities.remove_row(M, rrow=i)
                T = utilities.remove_row(T, rrow=i)
                break
    if verbose == 5:
        print("DBG: get_dc_guess(): M and T, after removing LD lines")
        print(M)
        print(T)

    # Remove empty columns:
    # If a column is empty, we have no guess regarding the corresponding
    # node. It makes the matrix singular. -> Remove the col & remember
    # that we are _not_ calculating a guess for it.
    removed_index = []
    for i in range(M.shape[1] - 1, -1, -1):
        if not M[:, i].any():
            M = utilities.remove_row_and_col(M, rrow=M.shape[0], rcol=i)
            removed_index.append(i)

    if verbose > 3:
        print("DBG: get_dc_guess(): M and T, after removing empty columns.")
        print(M)
        print("T\n", T)

    # Now, we have a set of equations to be solved.
    # There are three cases:
    # 1. The M matrix has a different number of rows and columns.
    #    We use the Moore-Penrose matrix inverse to get
    #    the shortest length least squares solution to the problem
    #          M*x + T = 0
    # 2. The matrix is square.
    #    It seems that if the circuit is not pathological,
    #    we are likely to find a solution (the matrix has det != 0).
    #    I'm not sure about this though.

    if M.shape[0] != M.shape[1]:
        Rp = np.dot(np.linalg.pinv(M), T)
    else:  # case M.shape[0] == M.shape[1], use normal
        if np.linalg.det(M) != 0:
            try:
                Rp = np.dot(np.linalg.inv(M), T)
            except np.linalg.linalg.LinAlgError:
                eig = np.linalg.eig(M)[0]
                cond = abs(eig).max() / abs(eig).min()
                if verbose:
                    print("cond=" + str(cond) + ". No guess.")
                return None
        else:
            if verbose:
                print("Guess matrix is singular. No guess.")
            return None

    # Now we want to:
    # 1. Add voltages for the nodes for which we have no clue to guess.
    # 2. Append to each vector of guesses the values for currents in
    #    voltage defined elem.
    # Both them are set to 0
    for index in removed_index:
        Rp = np.concatenate((np.concatenate((Rp[:index, 0].reshape((-1, 1)),
                                             np.zeros((1, 1))), axis=0),
                             Rp[index:, 0].reshape((-1, 1))), axis=0)
    # add the 0s for the currents due to the voltage defined
    # elements (we have no guess for those...)
    if v_eq > 0:
        Rp = np.concatenate((Rp, np.zeros((v_eq, 1))), axis=0)

    if verbose and verbose < 4:
        print("done.")
    if verbose > 3:
        print("Guess:")
        print(Rp)

    return Rp
