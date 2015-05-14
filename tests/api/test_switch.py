# -*- coding: iso-8859-1 -*-
# switch.py
# Implementation of the voltage controlled switch model
# Copyright 2013 Giuseppe Venturini
#
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

from __future__ import (unicode_literals, absolute_import,
                        division, print_function)

import pylab
import numpy as np
import numpy.random

import ahkab
from ahkab.switch import vswitch_model

cli = False

def test_switch():
    """Test switch API"""
    # This is a small test to check all is OK
    VT = 1.3 * np.random.randn()
    VH = abs(1 * np.random.randn())
    RON = abs(1e3 * np.random.randn())
    ROFF = abs(1e4 * np.random.randn())
    # VT = 0.; VH=1.; RON=100;
    if cli:
        print("Testing a switch with:")
        print("VT: %g\tVH: %g\tRON:%g\tROFF:%g" % (VT, VH, RON, ROFF))
    m = vswitch_model(name='test', VT=VT, VH=VH, RON=RON, ROFF=ROFF)
    VOs = [-5., -2.5, 0.0, 2.5, 5.]
    VIs = [-5., -2.5, 0.0, 2.5, 5.]
    VMAX = 5
    POINTS = 50

    class dev_class:
        pass
    device = dev_class()
    # TEST 1, VO
    vsweep = (2 * VMAX * np.arange(POINTS) / float(POINTS) - VMAX)
    vsweep = np.concatenate((vsweep, vsweep[::-1]))
    vsweep = vsweep.tolist()
    while True:
        if 0. in vsweep:
            i = vsweep.index(0.)
            vsweep = vsweep[:i] + vsweep[i+1:]
        else:
            break
    device.is_on = False
    if cli:
        pylab.hold(True)
    for vin in VIs:
        i = []
        go = []
        gon = []
        for vo in vsweep:
            i += [m.get_i((vo, vin), device)]
            go += [m.get_go((vo, vin), device)]
            gon += [i[-1] / (vo + 1e-12)]
        if cli:
            pylab.subplot(221)
            pylab.plot(vsweep, i, 'o-', label='VIN=%g' % (vin,), ms=3)
            pylab.ylabel('Output current [A]')
            pylab.xlabel('Output voltage [V]')
            pylab.subplot(222)
            pylab.plot(vsweep, go, 'o', color='g', ms=3, label='DATA VO=%g' % (vo,))
            pylab.plot(vsweep, gon, '-', label='TH VO=%g' % (vo,))
            pylab.ylabel('Output g [1/V]')
            pylab.xlabel('Output voltage [V]')
        assert np.allclose(go, gon, rtol=1e-3, atol=1e-7)
    for vo in VOs:
        i = []
        go = []
        gon = []
        for vin in vsweep:
            i += [m.get_i((vo, vin), device)]
            go += [m.get_go((vo, vin), device)]
            gon += [i[-1] / (vo + 1e-12)]
        if cli:
            pylab.subplot(223)
            pylab.plot(vsweep, i, 'o-', label='VO=%g' % (vo,), ms=3)
            pylab.xlabel('Input voltage [V]')
            pylab.ylabel('Output current [A]')
            pylab.subplot(224)
            pylab.plot(vsweep, go, 'o', color='g', label='DATA VO=%g' % (vo,), ms=3)
            pylab.plot(vsweep, gon, '-', label='TH VO=%g' % (vo,))
            pylab.ylabel('Output g [1/V]')
            pylab.xlabel('Input voltage [V]')
    if cli:
        for sbp in (221, 222, 223, 224):
            pylab.subplot(sbp)
            pylab.legend()

if __name__ == '__main__':
    cli = True
    test_switch()
    pylab.show()
