# -*- coding: iso-8859-1 -*-
# constants.py
# To hold the various constants required
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

"""Constants useful for building equations and expressions describing
semiconductor physics
"""


from __future__ import (unicode_literals, absolute_import,
                        division, print_function)
import math

#: The electron charge :math:`e`.
e = 1.60217646e-19
#: Simulation temperature in Kelvin degrees.
T = 300
#: Reference temperature in Kelvin degrees.
Tref = 300
#: The Boltzmann constant
k = 1.3806503e-23


def Vth(T=Tref):
    """The thermal voltage: :math:`kT/q`.

    **Parameters:**

    T : float, optional
        The temperature at which the thermal voltage is to be evaluated.
        If not specified, it defaults to ``constants.Tref``.

    **Returns:**

    vth : float
        The thermal voltage, :math:`kT/q`.
    """
    return k * T / e


class silicon:
    """Silicon class

    Access this class as ``constants.si``.

    **Attributes**

    **esi**: permittivity of silicon.

    **eox**: permittivity of silicon dioxide.

    """
    def __init__(self):
        self.esi = 104.5 * 10 ** -12  # F/m
        self.eox = 34.5 * 10 ** -12  # F/m

    def Eg(self, T=Tref):
        """Energy gap of silicon at temperature ``T``

        **Parameters:**

        T : float, optional
            The temperature at which the thermal voltage is to be evaluated.
            If not specified, it defaults to ``constants.Tref``.

        **Returns:**

        Eg : float
            The energy gap, expressed in electron-volt (eV).
        """

        return (1.16 - 0.000702 * T ** 2 / (1108 + T))  # eV

    def ni(self, T=Tref):
        """Intrinsic Silicon carrier concentration at temperature ``T``

        **Parameters:**

        T : float, optional
            The temperature at which the thermal voltage is to be evaluated.
            If not specified, it defaults to ``constants.Tref``.

        **Returns:**

        ni : float
            The intrinsec carrier concentration.
        """

        return 1.45 * 10 ** 16 * (T / Tref) * math.exp(self.Eg(Tref) / (2 * Vth(Tref)) - self.Eg(T) / (2 * Vth(T)))

#: Silicon class instantiated.
si = silicon()
