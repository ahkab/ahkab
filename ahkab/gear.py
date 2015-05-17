# -*- coding: iso-8859-1 -*-
# gear.py
# Gear Linear Multi-Step (LMS) DF
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

# This file should be conforming to a standard

""" About the method:
    This is an implicit method, it means that to compute dx(n+1)/dt the value of x in (n+1) is required.
    We don't know it, since it's our objective).
    This method, as all other implicit methods, allows us to write the derivative as:
     dx(n+1)/dt = x_coeff * x(n+1) + const                  (ii)

    The get_df method returns those two vectors.

    Gear's LMS interpolates the solution in a number of points equal to its order. Since it's a implicit
    method, one of these is x(n+1). The values x(n), x(n-1)... x(n-(order+2)) need to be supplied to the method.
    We can write x(t) as:
     x(t)   = a0 + a1*(t(n+1) -  t ) + a2*( t(n+1) - t )^2 + ...        (i)
    The equation has <order> a coefficients, which we need to determine.
    For this reason, we write a system of "order" equations in this way:
     x(n+1) = a0 + a1*(t(n+1) - t(n+1)) + a2*(t(n+1) - t(n+1))^2 + ...
     x(n)   = a0 + a1*(t(n+1) -  t(n) ) + a2*( t(n+1) - t(n) )^2 + ...
     x(n-1) = a0 + a1*(t(n+1) - t(n-1)) + a2*(t(n+1) - t(n-1))^2 + ...


    Which may be rewritten as:

      z = A * a

    z is the vector of known values of x
    A is a time dependant matrix
    a is a vector made of the a* coeffiecients

    We don't need to explicit ALL of the a* coeffiecients. What we are really looking for is the derivative of x
    in t(n+1), dx(n+1)/dt in short.
    If we differentiate the relation (i):
      dx(t)/dt = -a1 - 2*a2*( t(n+1) - t ) - 3*a3*( t(n+1) - t )^2 ...
    Which evaluated in t = t(n+1) gives:
      dx(n+1)/dt = -1 * a1

    Our objective is then a1.
    From the previsious system we write:
     a = A^-1 * z
    a1 is a[1,0], which may be extracted in this way:
     et = [0 1 0 0 0 0 ...] (order elements)
     a1 = et * a = et * A^-1 * z
    Because of the associative prperty of matrix multiplication, we can write:

     P  = et * A^-1
     a1 = P[1, :] * z

    But, we don't know z[0,0] = x(n+1), we can split the above relation:

     a1 = P[1, 0] * x(n+1) + P[1, 1:] * z[1:, 0]

    We arrived to the relation written above (ii)
     dx(n+1)/dt = x_coeff * x(n+1) + const = -1*a1

    So:
     x_coeff = -1 * P[1, 0]
     const   = -1 * P[1, 1:] * z[1:, 0]

    This module uses a faster way to compute the values that doesn't require to invert the matrix.
    Anyway, from a theorical point of view, the above applies.
      """

from __future__ import (unicode_literals, absolute_import,
                        division, print_function)

import numpy as np
import sys

from scipy.misc import factorial

from . import printing

order = None
# FAST = True


def is_implicit():
    return True


def has_ff():
    return True


def get_required_values():
    """This returns two python arrays built this way:
    [ max_order_of_x, max_order_of_dx ]
    Where:
    Both the values are int, or None
    if max_order_of_x is set to k, the df method needs all the x(n-i) values of x,
    where i<=k (the value the function assumed i+1 steps before the one we will ask for the derivative).
    The same applies to max_order_of_dx, but regards dx(n)/dt
    None means that NO value is required.

    The first array has to be used if no prediction is required, the second are the values needed for prediction.
    """
    # notice that: it returns the same values, it should be order-1 if no prediction is required BUT
    # we build the required values in a (hopefully) fast way, that requires
    # one point more.
    return ((order, None), (order, None))


def get_df(pv_array, suggested_step, predict=False):
    """The array must be built in this way:
    It has to be an array of arrays. Each of them has the following structure:

    [time, np_matrix, np_matrix]

    Hence the pv_array[k] element is made of:
    _ time is the time in which the solution is valid: t(n-k)
    _ The first np_matrix is x(n-k)
    _ The second is d(x(n-k))/dt
    Values that are not needed may be set to None and they will be disregarded.

    if predict == True, it needs one more point to give a prediction
    of x at the suggested step.

    Returns: None if the incorrect values were given, or quits.
    Otherwise returns an array:
    _ the [0] element is the np matrix of coeffiecients (Nx1) of x(n+1)
    _ the [1] element is the np matrix of constant terms (Nx1) of x(n+1)
    The derivative may be written as:
    d(x(n+1))/dt = ret[0]*x(n+1) + ret[1]"""

    if order is None:
        printing.print_general_error(
            "You must set Gear's order before using it! e.g. gear.order = 5")
        sys.exit(1)

    s = []
    s.append(0)
    for index in range(1, order + 2):
        s.append(suggested_step + pv_array[0][0] - pv_array[index - 1][0])

    # build e[k, i]
    e = np.zeros((order + 2, order + 2))
    for k_index in range(1, order + 2):
        for i_index in range(1, order + 2):
            if i_index == k_index:
                e[k_index, i_index] = 1
            else:
                e[k_index, i_index] = s[i_index] / (s[i_index] - s[k_index])

    alpha = np.zeros((1, order + 2))
    for k_index in range(1, order + 2):
        alpha[0, k_index] = 1.0
        for j_index in range(order + 1):
            alpha[0, k_index] = alpha[0, k_index] * e[k_index, j_index + 1]

    # build gamma
    gamma = np.zeros((1, order + 1))
    for k_index in range(1, order + 1):
        gamma[0, k_index] = alpha[0, k_index] * \
            ((1.0 / s[order + 1]) - (1.0 / s[k_index]))

    gamma[0, 0] = 0
    for index in range(1, order + 1):
        gamma[0, 0] = gamma[0, 0] - gamma[0, index]

    # values to be returned
    C1 = gamma[0, 0]

    C0 = np.zeros(pv_array[0][1].shape)
    for index in range(order):
        C0 = C0 + gamma[0, index + 1] * pv_array[index][1]

    x_lte_coeff = 0
    for k_index in range(1, order + 1):
        x_lte_coeff = x_lte_coeff + \
            (s[k_index] ** (order + 1)) * \
            (-1.0 * gamma[0, k_index] / gamma[0, 0])
    x_lte_coeff = ((-1.0) ** (order + 1)) * \
        (1.0/factorial(order + 1)) * x_lte_coeff

    if predict:
        predict_x = np.zeros(pv_array[0][1].shape)
        for index in range(1, order + 2):  # order
            predict_x = predict_x + alpha[0, index] * pv_array[index - 1][1]

        predict_lte_coeff = -1.0/factorial(order + 1)
        for index in range(1, order + 2):
            predict_lte_coeff = predict_lte_coeff * s[index]
        # print predict_lte_coeff
        # print x_lte_coeff
    else:
        predict_x = None
        predict_lte_coeff = None

    return C1, C0, x_lte_coeff, predict_x, predict_lte_coeff
