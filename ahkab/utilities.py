# -*- coding: iso-8859-1 -*-
# utilities.py
# Utilities file for Ahkab simulator
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

"""
This file holds miscellaneous utility functions needed by the simulator.
"""
from __future__ import unicode_literals, absolute_import

import os
import os.path
import operator

import numpy as np

from . import printing
from . import options

# this is the machine precision on my Intel x86
EPS = 2.22044604925e-16


def expand_matrix(matrix, add_a_row, add_a_col):
    """Adds a  row and/or a column to the given matrix.
    Args:
    matrix - the matrix to be manipulated
    add_a_row - boolean, if set to true adds a row
    add_a_col - boolean, if set to true adds a column

    Returns the new matrix"""
    (n_row, n_col) = matrix.shape
    if add_a_col:
        col = np.mat(np.zeros((n_row, 1)))
        matrix = np.concatenate((matrix, col), axis=1)
    if add_a_row:
        if add_a_col:
            n_col = n_col + 1
        row = np.mat(np.zeros((1, n_col)))
        matrix = np.concatenate((matrix, row), axis=0)
    return matrix


def remove_row_and_col(matrix, rrow=0, rcol=0):
    """Removes a row and a column from the matrix.
    rrow and rcol must be positive, or None is returned
    By default the first row and the first column are removed
    If you don't wish to remove one of them, supply a index that is greater
    than the matrix size.
    eg. matrix is 3x3, you want to remove just the second row of matrix, supply:
    rrow=1 and rcol=10 (or any number bigger than 2)
    """
    if rrow < 0 or rcol < 0:
        return_matrix = None
    else:
        return_matrix = np.vstack(
            (np.hstack((matrix[0:rrow, 0:rcol], matrix[0:rrow, rcol + 1:])), np.hstack((matrix[rrow + 1:, 0:rcol], matrix[rrow + 1:, rcol + 1:]))))
    return return_matrix


def remove_row(matrix, rrow=0):
    """Removes a row from a matrix.
    rrow is the index of the row to be removed.

    Returns: the matrix without the row, or none if rrow is invalid."""
    if rrow < 0 or rrow > matrix.shape[0] - 1:
        return_matrix = None
    else:
        return_matrix = np.vstack((matrix[:rrow, :], matrix[rrow + 1:, :]))
    return return_matrix


def check_file(filename):
    """Checks whether the supplied path refers to a valid file.
    Returns:
    True if it's found (and is a file)
    False, otherwise.
    """
    filename = os.path.abspath(filename)
    if not os.path.exists(filename):
        printing.print_general_error(filename + " not found.")
        ret = False
    elif not os.path.isfile(filename):
        printing.print_general_error(filename + " is not a file.")
        ret = False
    else:
        ret = True
    return ret

# Use scipy.factorial


def fact(num):
    """Returns: num!"""
    if num == 1:
        return 1
    return reduce(operator.mul, xrange(2, num + 1))


def calc_eps():
    """Returns the machine precision."""
    _eps = 1.0
    while(1 + _eps > 1):
        _eps = _eps / 2
    return _eps * 2


class combinations:

    """This class is an iterator that returns all the k-combinations
    _without_repetition_ of the elements of the supplied list.

    Each combination is made of a subset of the list, consisting of k
    elements.
    """

    def __init__(self, L, k):
        """This method initializes the class.
        L is the set from which the elements are taken.
        k is the size of the subset, the number of elements to be taken
        """
        self.L = L
        self.k = k
        self._sub_iter = None
        self._i = 0
        if len(self.L) < k:
            raise Exception, "The set has to be bigger than the subset."
        if k <= 0:
            raise Exception, "The size of the subset has to be positive."

    def __iter__(self):
        return self

    def next(self):
        """Get the next combination.
        Returns a list.
        """
        # It's recursive
        if self.k > 1:
            if self._sub_iter == None:
                self._sub_iter = combinations(self.L[self._i + 1:], self.k - 1)
            try:
                nxt = self._sub_iter.next()
                cur = self.L[self._i]
            except StopIteration:
                if self._i < len(self.L) - self.k:
                    self._i = self._i + 1
                    self._sub_iter = combinations(
                        self.L[self._i + 1:], self.k - 1)
                    return self.next()
                else:
                    raise StopIteration
        else:
            nxt = []
            if self._i < len(self.L):
                cur = self.L[self._i]
                self._i = self._i + 1
            else:
                raise StopIteration

        return [cur] + nxt


class log_axis_iterator:

    """This iterator provides the values for a logarithmic sweep.
    """

    def __init__(self, max, min, nsteps):
        self.inc = 10 ** ((np.log10(max) - np.log10(min)) / nsteps)
        self.max = max
        self.min = min
        self.index = 0
        self.current = min
        self.nsteps = nsteps

    def next(self):
        """Iterator method: get the next value
        """
        if self.index < self.nsteps:
            self.current = self.current * self.inc
            ret = self.current
        else:
            raise StopIteration
        self.index = self.index + 1
        return ret

    def __getitem__(self, i):
        """Iterator method: get a particular value (n. i)
        """
        if i == 0:
            ret = self.min
        elif i < self.nsteps:
            ret = self.min * self.inc ** i
        else:
            ret = None
        return ret

    def __iter__(self):
        """Required iterator method.
        """
        return self


class lin_axis_iterator:

    """This iterator provides the values for a linear sweep.
    """

    def __init__(self, max, min, nsteps):
        self.inc = (max - min) / nsteps
        self.max = max
        self.min = min
        self.index = 0
        self.current = min
        self.nsteps = nsteps

    def next(self):
        """Iterator method: get the next value
        """
        if self.index == 0:
            pass  # return min
        elif self.index < self.nsteps:
            self.current = self.current + self.inc
        else:
            raise StopIteration
        ret = self.current
        self.index = self.index + 1
        return ret

    def __getitem__(self, i):
        """Iterator method: get a particular value (n. i)
        """
        if i < self.nsteps:
            ret = self.min + self.inc * i
        else:
            ret = None
        return ret

    def __iter__(self):
        """Required iterator method.
        """
        return self


def Celsius2Kelvin(cel):
    return cel + 273.15


def Kelvin2Celsius(kel):
    return kel - 273.15

def convergence_check(x, dx, residuum, nv_minus_one, debug=False):
    if not hasattr(x, 'shape'):
        x = np.mat(np.array(x))
        dx = np.mat(np.array(dx))
        residuum = np.mat(np.array(residuum))
    vcheck, vresults = voltage_convergence_check(
        x[:nv_minus_one, 0], dx[:nv_minus_one, 0], residuum[:nv_minus_one, 0])
    icheck, iresults = current_convergence_check(
        x[nv_minus_one:], dx[nv_minus_one:], residuum[nv_minus_one:])
    return vcheck and icheck, vresults + iresults


def voltage_convergence_check(x, dx, residuum, debug=False):
    return custom_convergence_check(x, dx, residuum, er=options.ver, ea=options.vea, eresiduum=options.iea, debug=debug)


def current_convergence_check(x, dx, residuum, debug=False):
    return custom_convergence_check(x, dx, residuum, er=options.ier, ea=options.iea, eresiduum=options.vea, debug=debug)


def custom_convergence_check(x, dx, residuum, er, ea, eresiduum, vector_norm=lambda v: abs(v), debug=False):
    all_check_results = []
    if not hasattr(x, 'shape'):
        x = np.mat(np.array(x))
        dx = np.mat(np.array(dx))
        residuum = np.mat(np.array(residuum))
    if x.shape[0]:
        if not debug:
            ret = np.allclose(x, x + dx, rtol=er, atol=ea) and \
                np.allclose(residuum, np.zeros(
                               residuum.shape), atol=eresiduum, rtol=0)
        else:
            for i in range(x.shape[0]):
                if vector_norm(dx[i, 0]) < er * vector_norm(x[i, 0]) + ea and vector_norm(residuum[i, 0]) < eresiduum:
                    all_check_results.append(True)
                else:
                    all_check_results.append(False)
                if not all_check_results[-1]:
                    break

            ret = not (False in all_check_results)
    else:
        # We get here when there's no variable to be checked. This is because there aren't variables
        # of this type.
        # Eg. the circuit has no voltage sources nor voltage defined elements. In this case, the actual check is done
        # only by current_convergence_check, voltage_convergence_check always
        # returns True.
        ret = True

    return ret, all_check_results
