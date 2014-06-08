# -*- coding: iso-8859-1 -*-
# implicit_euler.py
# Implicit euler DF
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

""" This module implements IE (aka Backward Euler) and a first order prediction formula"""

import numpy as np

order = 1


def is_implicit():
    return True


def get_df_coeff(step):
    return (1.0 / step, -1.0 / step)


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
    # we need just the x(n) value
    return ((0, None), (1, None))


def has_ff():
    return True
# def get_df_lte_coeff(pv_array, suggested_step):
    # return 1 / 2 * suggested_step


def get_df(pv_array, suggested_step, predict=True):
    """The array must be built in this way:
    It must be an array of these array:

    [time, np_matrix, np_matrix]

    Hence the pv_array[k] element is made of:
    _ time is the time in which the solution is valid: t(n-k)
    _ The first np_matrix is x(n-k)
    _ The second is d(x(n-k))/dt
    Values that are not needed may be None, they will be disregarded
    Returns None if the incorrect values were given.
    Otherwise returns an array:
    _ the [0] element is the coeffiecient of x(n+1) (scalar)
    _ the [1] element is the np matrix of constant terms (Nx1) of x(n+1)
    _ the [2] element is the coefficient of the lte of x(n+1) (scalar)
    _ the [2] element is the predicted value of x(n+1) (np matrix), if predicted == True
      otherwise it's None
    _ the [3] element is the coefficient of the lte of the prediction (np matrix), if predict,
      otherwise None
    The derivative may be written as:
    d(x(n+1))/dt = ret[0]*x(n+1) + ret[1]"""

    # print pv_array
    # we'll use just the first column, since our method needs just x(n)

    if len(pv_array[0]) != 3:
        return None
    if pv_array[0][1] is None:
        return None

    # xold = np.mat(pv_array[0][1].copy())
    x_lte_coeff = 0.5 * suggested_step

    if predict and len(pv_array) > 1 and pv_array[1][1] is not None:
        predict = np.mat(np.zeros(pv_array[0][1].shape))
        for index in xrange(predict.shape[0]):
            predict[index, 0] = (pv_array[0][1][index, 0] - pv_array[1][1][index, 0]) \
                / (pv_array[0][0] - pv_array[1][0]) * suggested_step + pv_array[0][1][index, 0]
        predict_lte_coeff = -0.5 * suggested_step * \
            (pv_array[0][0] + suggested_step - pv_array[1][0])
    else:
        predict = None
        predict_lte_coeff = None

    # print x_lte_coeff
    # print predict_lte_coeff

    return (1.0 / suggested_step, -1.0 / suggested_step * pv_array[0][1], x_lte_coeff, predict, predict_lte_coeff)
