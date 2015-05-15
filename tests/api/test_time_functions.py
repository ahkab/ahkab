# -*- coding: iso-8859-1 -*-
# test_time_functions.py
# Unit tests for the time functions
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
# along with ahkab. If not, see <http://www.gnu.org/licenses/>.

from __future__ import unicode_literals, print_function, division

import numpy as np
import sympy

from sympy.functions.special.delta_functions import Heaviside

from ahkab import time_functions

def test_sffm():
    """Test time_functions.sffm"""
    vo, va, fc, mdi, fs, td, time = sympy.symbols('vo, va, fc, mdi, fs, td,' +
                                                 ' time')
    F = vo + va*sympy.sin(2*sympy.pi*fc*(time - td) +
                         mdi*sympy.sin(2*sympy.pi*fs*(time - td)))* \
            Heaviside(time - td)
    von, van, fcn, mdin, fsn, tdn = 0.1, 1e-3, 20e3, 5, 1e3, 1e-3
    f = time_functions.sffm(vo=von, va=van, fc=fcn, mdi=mdin, fs=fsn, td=tdn)
    FS = sympy.lambdify(time, sympy.N(F.subs(dict(vo=von, va=van, fc=fcn,
                                                 mdi=mdin, fs=fsn,
                                                 td=tdn))))
    t = np.linspace(0, 1e-3, 3e3)
    for ti in t:
        assert np.allclose(f.value(ti), float(FS(ti)), rtol=1e-4)

def test_am():
    """Test time_functions.am"""
    sa, fc, fm, oc, td, time = sympy.symbols('sa, fc, fm, oc, td, time')
    F = sa*(oc + sympy.sin(2*sympy.pi*fm*(time - td)))* \
        sympy.sin(2*sympy.pi*fc*(time - td)) * \
        Heaviside(time - td)

    san, fcn, fmn, ocn, tdn = 10., 1., 1e3, 100, 1e-3
    f = time_functions.am(sa=san, fc=fcn, fm=fmn, oc=ocn, td=tdn)
    FS = sympy.lambdify(time, sympy.N(F.subs(dict(sa=san, fc=fcn, fm=fmn,
                                                  oc=ocn, td=tdn))))
    t = np.linspace(0, 20e-3, 3e3)
    for ti in t:
        assert np.allclose(f.value(ti), float(FS(ti)), rtol=1e-4)


if __name__ == '__main__':
    test_sffm()
