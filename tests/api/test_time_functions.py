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
        assert np.allclose(f(ti), float(FS(ti)), rtol=1e-4)

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
        assert np.allclose(f(ti), float(FS(ti)), rtol=1e-4)

def test_sin():
    """Test time_functions.sin"""
    vo, va, freq, td, theta, phi, time = sympy.symbols('vo, va, freq, td, ' +
                                                       'theta, phi, time')
    F = vo + va*sympy.sin(sympy.pi*phi/180)*Heaviside(td - time) + \
             va*sympy.sin(2*sympy.pi*freq*(time - td) + sympy.pi*phi/180)* \
             sympy.exp(-(time - td)*theta)*Heaviside(time - td)

    von, van, freqn, tdn, thetan, phin = 0.1, 3, 20e3, 0.5e-3, 1e2, 50.
    f = time_functions.sin(vo=von, va=van, freq=freqn, td=tdn, theta=thetan,
                           phi=phin)
    FS = sympy.lambdify(time, sympy.N(F.subs(dict(vo=von, va=van, freq=freqn,
                                                  td=tdn, theta=thetan,
                                                  phi=phin))))
    t = np.linspace(0, 1e-3, 3e3)
    for ti in t:
        assert np.allclose(f(ti), float(FS(ti)), rtol=1e-4)

def test_exp():
    """Test time_functions.exp"""
    v1, v2, td1, tau1, td2, tau2, time = sympy.symbols('v1, v2, td1, tau1, ' +
                                                       'td2, tau2, time')
    F = v1 + (v2 - v1)*(1 - sympy.exp(-(time - td1)/tau1))*Heaviside(time - td1) +\
             (v1 - v2)*(1 - sympy.exp(-(time - td2)/tau2))*Heaviside(time - td2)
    v1n, v2n, td1n, tau1n, td2n, tau2n = 0.1, 3, 0.33e-3, 0.2e-3, .66e-3, .11e-3
    f = time_functions.exp(v1=v1n, v2=v2n, td1=td1n, tau1=tau1n, td2=td2n, tau2=tau2n)
    FS = sympy.lambdify(time, sympy.N(F.subs(dict(v1=v1n, v2=v2n, td1=td1n,
                                                  tau1=tau1n, td2=td2n,
                                                  tau2=tau2n))))
    t = np.linspace(0, 1e-3, 3e3)
    for ti in t:
        assert np.allclose(f(ti), float(FS(ti)), rtol=1e-4)

def test_pulse():
    """Test time_functions.pulse"""
    # as we cannot define a periodic function in sympy, we will only check one
    # period
    v1, v2, td, tr, pw, tf, per, time = sympy.symbols('v1, v2, td, tr, pw, ' +
                                                      'tf, per, time')
    F = v1*Heaviside(td - time) + \
        ((v2 - v1)/tr*time + v1 +(v1 - v2)*td/tr)* \
        Heaviside(time - td)*Heaviside(td + tr - time) + \
        v2*Heaviside(time - td - tr)*Heaviside(td + tr + pw - time) + \
        ((v1 - v2)/tf*time + v2 + (-v1*(pw + td + tr) + v2*(pw + td + tr))/tf)* \
        Heaviside(time - td - tr - pw)*Heaviside(td + tr + pw + tf - time) + \
        v1*Heaviside(time - td - tr - pw - tf)

    v1n, v2n, tdn, trn, pwn, tfn, pern = (-2, 3, 0.1e-3, 0.05e-3, .5e-3, .01e-3,
                                         1e-3)
    f = time_functions.pulse(v1=v1n, v2=v2n, td=tdn, tr=trn, pw=pwn, tf=tfn,
                             per=pern)
    FS = sympy.lambdify(time, sympy.N(F.subs(dict(v1=v1n, v2=v2n, td=tdn,
                                                  tr=trn, pw=pwn, tf=tfn,
                                                  per=pern))))
    t = np.linspace(0, 1e-3, 1e3)
    for ti in t:
        assert np.allclose(f(ti), float(FS(ti)), rtol=1e-4)

def test_pwl():
    """Test time_functions.pwl"""
    # we define a pulse shape, with repeat and we check that it matches with
    # its sympy implementation
    # as we cannot define a periodic function in sympy, we will use a small hack
    v1, v2, td, tr, pw, tf, per, time = sympy.symbols('v1, v2, td, tr, pw, ' +
                                                      'tf, per, time')
    F = v1*Heaviside(td - time) + \
        ((v2 - v1)/tr*time + v1 +(v1 - v2)*td/tr)* \
        Heaviside(time - td)*Heaviside(td + tr - time) + \
        v2*Heaviside(time - td - tr)*Heaviside(td + tr + pw - time) + \
        ((v1 - v2)/tf*time + v2 + (-v1*(pw + td + tr) + v2*(pw + td + tr))/tf)* \
        Heaviside(time - td - tr - pw)*Heaviside(td + tr + pw + tf - time) + \
        v1*Heaviside(time - td - tr - pw - tf)

    v1n, v2n, tdn, trn, pwn, tfn, pern = (-2, 3, 0.1e-3, 0.05e-3, .5e-3, .01e-3,
                                         1e-3)
    x = [0., 0.05e-3, 0.55e-3, 0.56e-3, 1e-3]
    y = [-2, 3, 3, -2, -2]
    f = time_functions.pwl(x, y, repeat=True, td=0.1e-3, repeat_time=1e-3)
    FS = sympy.lambdify(time, sympy.N(F.subs(dict(v1=v1n, v2=v2n, td=tdn,
                                                  tr=trn, pw=pwn, tf=tfn,
                                                  per=pern))))
    t = np.linspace(0, 1e-3, 1e3)
    for ti in t:
        assert np.allclose(f(ti), float(FS(ti)), rtol=1e-4)
    for ti in t:
        assert np.allclose(f(ti+1e-3), float(FS(ti)), rtol=1e-4)

if __name__ == '__main__':
    test_sffm()
    test_am()
    test_sin()
    test_exp()
    test_pulse()
    test_pwl()
