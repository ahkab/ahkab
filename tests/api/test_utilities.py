# -*- coding: iso-8859-1 -*-
# test_utilities.py
# Utilities module test file for the Ahkab simulator
# Copyright 2015 Giuseppe Venturini

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
This module contains test functions for the utility module.

"""

from __future__ import (unicode_literals, absolute_import,
                        division, print_function)

import os, os.path
from nose.tools import raises

import ahkab
from ahkab import utilities
from ahkab.utilities import *

import numpy as np

#: The machine epsilon, the upper bound on the relative error due to rounding in
#: floating point arithmetic.
def test_eps():
    """Test utilities.EPS"""
    assert utilities.EPS == np.finfo(float).eps

def test_expand_matrix():
    """Test utilities.expand_matrix()"""
    m = np.eye(5)
    assert not (m - utilities.expand_matrix(m)).any()
    assert (6, 5) == utilities.expand_matrix(m, add_a_row=True).shape
    assert (5, 6) == utilities.expand_matrix(m, add_a_col=True).shape
    assert (6, 6) == utilities.expand_matrix(m, add_a_col=True, add_a_row=True).shape

def test_set_submatrix():
    """Test utilities.set_submatrix()"""
    m1 = np.array(((1, 1, 1), (2, 2, 2), (3, 3, 3)))
    m2 = np.eye(10)
    set_submatrix(row=3, col=3, dest_matrix=m2, source_matrix=m1)
    assert not (m2 - \
    np.array([[ 1.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
           [ 0.,  1.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
           [ 0.,  0.,  1.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
           [ 0.,  0.,  0.,  1.,  1.,  1.,  0.,  0.,  0.,  0.],
           [ 0.,  0.,  0.,  2.,  2.,  2.,  0.,  0.,  0.,  0.],
           [ 0.,  0.,  0.,  3.,  3.,  3.,  0.,  0.,  0.,  0.],
           [ 0.,  0.,  0.,  0.,  0.,  0.,  1.,  0.,  0.,  0.],
           [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  1.,  0.,  0.],
           [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  1.,  0.],
           [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  1.]])).any()

def test_remove_row_and_col():
    """Test utilities.remove_row_and_col()"""
    m1 = np.eye(5)
    assert not (remove_row_and_col(m1, rrow=0, rcol=0) - np.eye(4)).any()
    assert remove_row_and_col(m1, rrow=None, rcol=0).shape == (5,4)
    assert remove_row_and_col(m1, rrow=0, rcol=None).shape == (4,5)

def test_remove_row():
    """Test utilities.remove_col()"""
    m1 = np.eye(5)
    assert not (remove_row(m1, rrow=0) - remove_row_and_col(m1, rrow=0, rcol=None)).any()
    assert remove_row(m1, rrow=0).shape == (4,5)

def test_check_file1():
    """Test utilities.check_file() 1/3"""
    usf = os.path.realpath(__file__)
    assert check_file(usf)

@raises(IOError)
def test_check_file2():
    """Test utilities.check_file() 2/3"""
    usf = os.path.realpath(__file__)
    assert not check_file(usf+'12345')

@raises(IOError)
def test_check_file3():
    """Test utilities.check_file() 3/3"""
    usf = os.path.realpath(__file__)
    usd = os.path.dirname(usf)
    assert not check_file(usd)

def test_tuplinator():
    """Test utilities.tuplinator()"""
    a = [[[1,2,3], [1, []]]]
    assert tuplinator(a) == (((1, 2, 3), (1, ())),)

def test_combinations():
    """Test utilities.combinations"""
    a = combinations([1, 2, 3, 4, 5, 6], 2)
    a = list(a)
    d = {(1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (2, 3), (2, 4), (2, 5), (2, 6),
         (3, 4), (3, 5), (3, 6), (4, 5), (4, 6), (5, 6)}
    assert d == set(tuplinator(a))

def test_log_axis_iterator():
    """Test utilities.log_axis_iterator"""
    # logspace with endpoint
    a = np.logspace(0, 3, 1000, True)
    # iterator to list to array
    b = np.array(list(log_axis_iterator(1e3, 1, 1000)))
    assert abs((a - b).mean()) < .2


def test_lin_axis_iterator():
    """Test utilities.lin_axis_iterator()"""
    # logspace with endpoint
    a = np.linspace(1, 3, 1000, True)
    # iterator to list to array
    b = np.array(list(lin_axis_iterator(3, 1, 1000)))
    assert abs((a - b).mean()) < 2e-3


def test_Celsius2Kelvin():
    """Tests utilities.Celsius2Kelvin()"""
    assert Celsius2Kelvin(-273.15) == 0


def test_Kelvin2Celsius():
    """Test utilities.Kelvin2Celsius()"""
    assert Kelvin2Celsius(273.15) == 0

def test_custom_convergence_check():
    """Test utilities.custom_array_convergence_check()"""
    x = np.array(((1, 1, 1, 1),), dtype=np.float64)
    dx = np.array(((.1, .01, .001, .0001),), dtype=np.float64)
    residuum = np.array(((1e-6, 1e-6, 1e-6, 1e-6),), dtype=np.float64)
    # pass
    assert custom_convergence_check(x, dx, residuum, er=1, ea=.2, eresiduum=1e-5, debug=False)[0]
    assert not custom_convergence_check(x, dx, residuum, er=.01, ea=0, eresiduum=1e-5, debug=False)[0]
    assert not custom_convergence_check(x, dx, residuum, er=.01, ea=.02, eresiduum=1e-5, debug=False)[0]
    assert not custom_convergence_check(x, dx, residuum, er=.01, ea=.2, eresiduum=1e-7, debug=False)[0]

    assert custom_convergence_check(x, dx, residuum, er=1, ea=.2, eresiduum=1e-5, debug=True)[0]
    assert not custom_convergence_check(x, dx, residuum, er=.01, ea=0, eresiduum=1e-5, debug=True)[0]
    assert not custom_convergence_check(x, dx, residuum, er=.01, ea=.02, eresiduum=1e-5, debug=True)[0]
    assert not custom_convergence_check(x, dx, residuum, er=.01, ea=.2, eresiduum=1e-7, debug=True)[0]

def test_check_step_and_points1():
    """Test check_step_and_points() 1/4"""
    a = check_step_and_points(step=1., points=100, period=100.)
    assert a == (101, 1.)

def test_check_step_and_points2():
    """Test check_step_and_points() 2/4"""
    a = check_step_and_points(step=1., points=None, period=100.)
    assert a == (101, 1.)

def test_check_step_and_points3():
    """Test check_step_and_points() 3/4"""
    a = check_step_and_points(step=None, points=101, period=100.)
    assert a == (101, 1.)

def test_check_step_and_points4():
    """Test check_step_and_points() 4/4"""
    a = check_step_and_points(step=1., points=None, period=100.,
                              default_points=101)
    assert a == (101, 1.)

def test_check_circuit1():
    """Test utilities.check_circuit() 1/4"""
    c = ahkab.Circuit('New')
    c.add_resistor('R1', '0', '0', 1e3)
    v, _ = ahkab.utilities.check_circuit(c)
    assert not v

def test_check_circuit2():
    """Test utilities.check_circuit() 2/4"""
    c = ahkab.Circuit('New')
    c.add_resistor('R1', 'n1', 'n2', 1e3)
    c.add_resistor('R2', 'n1', 'n2', 1e3)
    v, _ = ahkab.utilities.check_circuit(c)
    assert not v

def test_check_circuit3():
    """Test utilities.check_circuit() 3/4"""
    c = ahkab.Circuit('New')
    c.add_resistor('R1', 'n1', 'gnd', 1e3)
    v, _ = ahkab.utilities.check_circuit(c)
    assert not v

def test_check_circuit4():
    """Test utilities.check_circuit() 4/4"""
    c = ahkab.Circuit('New')
    c.add_resistor('R1', 'n1', 'n2', 1e3)
    c.add_resistor('R1', 'n2', '0', 1e3)
    v, _ = ahkab.utilities.check_circuit(c)
    assert not v

def test_memoize():
    """Test utilities.memoize()"""
    @ahkab.utilities.memoize
    def test(a):
        return a**2
    # we just check that cache and misses
    # return the same values
    assert test(5) == test(5)

