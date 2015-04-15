# -*- coding: iso-8859-1 -*-
# trap.py
# Trap DF (with a second order prediction)
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
This file implements the TRAP DF and a second order prediction formula.
"""

from __future__ import (unicode_literals, absolute_import,
                        division, print_function)

import numpy as np
from numpy.linalg import inv

order = 2


def is_implicit():
    """Returns: boolean"""
    return True


def get_required_values():
    """This returns a python array built this way:
    [ max_order_of_x, max_order_of_dx ]
    Where:
    Both the values are int, or None
    if max_order_of_x is set to k, the df method needs all the x(n-i) values of x,
    where i<=k (the value the function assumed i+1 steps before the one we will ask for the derivative).
    The same applies to max_order_of_dx, but regards dx(n)/dt
    None means that NO value is required.

    The first array has to be used if no prediction is required, the second are the values needed for prediction.
    """
    # we need x(n) and dx(n)/dt
    return ((0, 0), (2, None))


def has_ff():
    """Has the method a Forward Formula for prediction?
    Returns: boolean
    """
    return True


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
    _ the [0] element is the np matrix of coeffiecients (Nx1) of x(n+1)
    _ the [1] element is the np matrix of constant terms (Nx1) of x(n+1)
    The derivative may be written as:
    d(x(n+1))/dt = ret[0]*x(n+1) + ret[1]"""

    # our method needs x(n) dx(n)/dt
    if len(pv_array[0]) != 3:
        return None
    if pv_array[0][1] is None or pv_array[0][2] is None:
        return None

    C1 = 2.0 / suggested_step
    C0 = -1.0 * ((2.0 / suggested_step) * pv_array[0][1] + pv_array[0][2])
    x_lte_coeff = (1.0 / 12) * (suggested_step ** 3)

    if predict and len(pv_array) > 2 and pv_array[0][1] is not None and pv_array[1][1] is not None and \
            pv_array[2][1] is not None:
        A = np.mat(np.zeros((len(pv_array[0]), len(pv_array[0]))))
        A[0, :] = 0
        A[:, 0] = 1
        for row in range(1, A.shape[0]):
            for col in range(1, A.shape[0]):
                A[row, col] = (pv_array[row][0] - pv_array[0][0]) ** (col)
        Ainv = inv(A)

        predict_x = np.mat(np.zeros(pv_array[0][1].shape))
        z = np.mat(np.zeros((3, 1)))
        for var in range(pv_array[0][1].shape[0]):

            for index in range(z.shape[0]):
                z[index, 0] = pv_array[index][1][var, 0]
            alpha = Ainv * z
            predict_x[var, 0] = alpha[2, 0] * (suggested_step ** 2) + \
                                alpha[1, 0] * suggested_step + alpha[0, 0]
        predict_lte_coeff = (-1.0 / 6.0) * suggested_step * \
                            (pv_array[0][0] + suggested_step - pv_array[1][0]) * \
                            (pv_array[0][0] + suggested_step - pv_array[2][0])
    else:
        predict_x, predict_lte_coeff = (None, None)
    return [C1, C0, x_lte_coeff, predict_x, predict_lte_coeff]
