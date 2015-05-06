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

"""This module implements the Implicit Euler (IE, aka Backward Euler, BE) and a
first-order forward formula (FF) to be used for prediction.

The formula is:

.. math::

    x'_{n+1} = C_0 x_{n+1} + C_1 x_{n}

Where:

* :math:`C_0 = 1/h`
* :math:`C_1 = -1/h`

The backward Euler method is not only A-stable, making it suitable for the solution of
stiff equations but is even L-stable.

Module reference
================

"""

from __future__ import (unicode_literals, absolute_import,
                        division, print_function)

import numpy as np

#: The order of the differentiation formula
order = 1

def is_implicit():
    """Is this differentiation formula implicit?"""
    return True


def get_df_coeff(step):
    """Get the coefficients for a Backward Euler differentiation step

    The first coefficient is the factor for the new point :math:`x_{n+1}`,
    the second is the one for the previous point :math:`x_{n}`.

    If the step value is :math:`h`, this method returns:

    .. math::

        [1/h, -1/h]

    **Parameters:**

    step : float
        The differentiation formula step value.

    **Returns:**

    c0, c1 : floats
        The coefficients of :math:`x_{n+1}` and :math:`x_{n}`.
    """
    return 1.0/step, -1.0/step


def get_required_values():
    """Get what values are required by the DF and the FF

    **Returns**

    The method returns two tuples, each of them having the form:

    ``[max_order_of_x, max_order_of_dx]``

    The first tuple is the one to be considered if no Forward Formula (FF) is
    needed, the second if the FF is also required.

    Both the values in each tuple can be either of type ``int`` or be set to
    ``None``.

    If ``max_order_of_x`` is set to an arbitrary positive integer value
    :math:`k`, the Differentiation Formula (DF) needs all the :math:`x_{n-i}`
    values of :math:`x`, where :math:`i \le k` (the value x has :math:`i+1`
    steps before the one we will ask for the derivative).  The same applies to
    ``max_order_of_dx``, but it regards :math:`dx/dt` instead of :math:`x`.

    If ``max_order_of_x`` or ``max_order_of_dx`` are set to ``None``, that means
    that no value of :math:`x`, or :math:`dx/dt`, is required.

    In the case at hand, where the formula is the Backward Euler (BE, aka
    Implicit Euler, IE), this method will return:

    ``((0, None), (1, None))``

    """
    # we need just the x(n) value
    return ((0, None), (1, None))


def has_ff():
    """Is a forward formula for prediction available?

    **Returns:**

    ``True``

    """
    return True


def get_df(pv_array, suggested_step, predict=True):
    """Get the coefficients for the DF, FF and LTE calculation

    **Parameters:**

    pv_array : list
        It must be an list of lists, each of them having the structure
        ``[time, xnk, dxnk]``.

        In particular, the ``pv_array[k]`` element of ``pv_array`` is composed
        of:

        * ``time``, float, which is the time at which the solution is valid:
          ``t(n-k)``,
        * ``xnk``, ndarray, which is :math:`x_{n-k}`,
        * ``dxnk``, ndarray, :math:`dx_{n-k}/dt`.

        The length of ``pv_array`` has to match the value returned by
        :func:`get_required_values`.

        Any values that are not needed may be set to ``None``, and they will be
        disregarded.

    suggested_step : float
        The step that is (expected) to be used in the DF. It is only an
        expectation because it may be rejected at a later stage if there is step
        control enabled.
    predict : boolean, optional
        Whether the terms for a prediction formula are required as well or not.
        Defaults to ``True``.

    **Returns:**

    ret : tuple
        ``ret`` is a tuple of 5 elements, where:

        * the [0] element is the coeffiecient of :math:`x_{n+1}` (scalar),
        * the [1] element is the matrix of constant terms of shape ``(Nx1)`` of
          :math:`x_{n+1}`,
        * the [2] element is the coefficient of the LTE of :math:`x_{n+1}` (scalar),
        * the [3] element is the predicted value of :math:`x_{n+1}` (matrix), only
          available if the ``predict`` parameter is set to ``True``. Otherwise it's
          ``None``.
        * the [4] element is the coefficient of the LTE of the prediction (matrix),
          also only available if the ``predict`` parameter is set to ``True``,
          otherwise, it is ``None``.

    .. note::

        With the returned values, the derivative may then be written as:

        .. math::

            \\frac{dx_{n+1}}{dt} = \\mathrm{ret[0]}\; x_{n+1} + \\mathrm{ret[1]}

    """

    # we'll use just the first column, since our method needs just x(n)
    if len(pv_array[0]) != 3:
        raise ValueError('Malformed pv_array: length == %d (3)' % len(pv_array))
    if pv_array[0][1] is None:
        raise ValueError('Malformed pv_array: pv_array[0][1] is None')

    x_lte_coeff = 0.5*suggested_step

    if predict and len(pv_array) > 1 and pv_array[1][1] is not None:
        predict = np.zeros(pv_array[0][1].shape)
        for index in range(predict.shape[0]):
            predict[index, 0] = (pv_array[0][1][index, 0] -
                                 pv_array[1][1][index, 0]) \
                                / (pv_array[0][0] -
                                   pv_array[1][0])*suggested_step \
                                + pv_array[0][1][index, 0]
        predict_lte_coeff = -0.5*suggested_step*\
                            (pv_array[0][0] + suggested_step - pv_array[1][0])
    else:
        predict = None
        predict_lte_coeff = None

    return (1./suggested_step, -1./suggested_step*pv_array[0][1],
            x_lte_coeff, predict, predict_lte_coeff)
