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
This file implements the Trapezoidal (TRAP) Differentiation Formula (DF) and a
second order prediction formula.

Module reference
----------------
"""

from __future__ import (unicode_literals, absolute_import,
                        division, print_function)

import numpy as np
from numpy.linalg import inv

from .py3compat import range_type

order = 2


def is_implicit():
    """Is this Differentiation Formula (DF) implicit?

    **Returns:**

    isit : boolean
        In this case, that's ``True``.
    """
    return True


def get_required_values():
    """Get info regarding what values are needed by the DF

    **Returns:**

    tpl : tuple of tuples
        A tuple of two tuples.

        * The first tuple indicates what past values of the unknown are needed
          for the DF.
        * The second tuple indicates what past values of the unknown are needed
          for the prediction method.

        In particular, each of the sub-tuples is built this way:
        ::

            (max_order_of_x, max_order_of_dx)

        Where both the values are either ``int``, or ``None``.
        If ``max_order_of_x`` is set to an integer value :math:`k`, the DF needs
        all the :math:`x_{n-i}` values of x, for all :math:`0 \\le i \\le k`. In
        the previous text, :math:`x_{n-i}` is the value the :math:`x` array
        assumed :math:`i` steps before the one we are considering for the
        derivative.

        Similar considerations apply to ``max_order_of_dx``, but regard rather
        :math:`dx_n/dt` instead of :math:`x_n`.

        If any of the values is set to ``None``, it is to be assume that no
        value is required.

    The first array has to be used if no prediction is required, the second are the values needed for prediction.
    """
    # we need x(n) and dx(n)/dt
    return ((0, 0), (2, None))


def has_ff():
    """Has the method a Forward Formula for prediction?

    **Returns:**

    doesit : bool
        In this particular case, this function always returns ``True``.
    """
    return True


def get_df(pv_array, suggested_step, predict=True):
    """Get the coefficients for DF and prediction formula

    **Parameters:**

    pv_array : sequence of sequences
        Each element of ``pv_array`` must be of the form:

        ::

            (time, x, derivate(x))

        In particular, the :math:`k` element of `pv_array` contains the values
        of:

        * :math:`t_{n-k}` (the time),
        * :math:`x_{n-k}`,
        * :math:`dx_{n-k}/dt`

        evaluated :math:`k` time steps before the current one, labeled
        :math:`n+1`.

        How many samples are necessary is given by
        :func:`ahkab.trap.get_required_values`.

        Values that are not needed may be set to ``None``, as they will be
        disregarded.
    suggested_step : float
        The step that will be used for the current iteration, provided the error
        will be deemed acceptable.
    predict : bool, optional
        Whether a prediction for :math:`x_n` is needed as well or not. Defaults
        to ``True``.

    **Returns:**

    ret : tuple
        The return value has the form:

        ::

            (C1, C0, x_lte_coeff, predict_x, predict_lte_coeff)

        The derivative may be written as:

        .. math::

            d(x(n+1))/dt = C1 x(n+1) + C0

        `x_lte_coeff` is the coefficient of the Local Truncation Error,
        `predict_x` is the predicted value for :math:`x` and `predict_lte_coeff`
        is the LTE coefficient for the prediction.

    :raises ValueError: if the `pv_array` is malformed.
    """

    # our method needs x(n) dx(n)/dt
    if len(pv_array[0]) != 3:
        raise ValueError('trap.get_df got a pv_array of wrong dimensions')
    if pv_array[0][1] is None or pv_array[0][2] is None:
        raise ValueError('trap.get_df got a pv_array that does not contain'
                         ' required values')
    if predict and (len(pv_array) < 2 or pv_array[0][1] is None or
        pv_array[1][1] is None or pv_array[2][1] is None):
        raise ValueError('trap.get_df got predict=True but pv_array that does'
                         ' not contain the required values')

    C1 = 2.0 / suggested_step
    C0 = -1.0 * (2.0 / suggested_step * pv_array[0][1] + pv_array[0][2])
    x_lte_coeff = 1.0 / 12 * suggested_step ** 3

    if predict and len(pv_array) > 2 and pv_array[0][1] is not None and \
        pv_array[1][1] is not None and pv_array[2][1] is not None:
        A = np.zeros((len(pv_array[0]), len(pv_array[0])))
        A[0, :] = 0
        A[:, 0] = 1
        for row in range_type(1, A.shape[0]):
            for col in range_type(1, A.shape[0]):
                A[row, col] = (pv_array[row][0] - pv_array[0][0]) ** col
        Ainv = inv(A)

        predict_x = np.zeros(pv_array[0][1].shape)
        z = np.zeros((3, 1))
        for var in range_type(pv_array[0][1].shape[0]):
            for index in range(z.shape[0]):
                z[index, 0] = pv_array[index][1][var, 0]
            alpha = np.dot(Ainv, z)
            predict_x[var, 0] = alpha[2, 0] * suggested_step**2 + \
                                alpha[1, 0] * suggested_step + alpha[0, 0]
        predict_lte_coeff = -1.0 / 6.0 * suggested_step * \
                            (pv_array[0][0] + suggested_step - pv_array[1][0]) * \
                            (pv_array[0][0] + suggested_step - pv_array[2][0])
    else:
        predict_x, predict_lte_coeff = (None, None)
    return C1, C0, x_lte_coeff, predict_x, predict_lte_coeff

