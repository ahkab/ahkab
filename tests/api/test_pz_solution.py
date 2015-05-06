# -*- coding: iso-8859-1 -*-
# test_pz_solution.py
# Unit tests for the PZ solution class
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

from __future__ import print_function, division, unicode_literals

import numpy as np
import ahkab

from nose import SkipTest
from ahkab import py3compat

def test_pz_solution():
    """Test results.pz_solution"""
    # Numeric test, not to run on PYPY
    if py3compat.PYPY:
        raise SkipTest
    bpf = ahkab.Circuit('RLC bandpass')
    bpf.add_inductor('L1', 'in', 'n1', 1e-6)
    bpf.add_capacitor('C1', 'n1', 'out', 2.2e-12)
    bpf.add_resistor('R1', 'out', bpf.gnd, 13)
    bpf.add_vsource('V1', 'in', bpf.gnd, dc_value=1, ac_value=1)

    pza = ahkab.new_pz('V1', ('out', bpf.gnd), x0=None, shift=1e3)
    r = ahkab.run(bpf, pza)['pz']

    str(r)

    np.allclose(r['p0']+r['p1'], -1034507*2, rtol=1e-3)
    np.allclose(abs(r['p0']-r['p1']), 107297253*2, rtol=1e-3)
    np.allclose(r['z0'], 0, rtol=1.)

    assert set(r.keys()) == {u'p0', u'p1', u'z0'}
    assert r['p0'] == r.get('p0')

    for i in {u'p0', u'p1', u'z0'}:
        assert i in r
        assert r.has_key(i)

    set(list(zip(*r.items()))[0]) == {u'p0', u'p1', u'z0'}
    set(list(zip(*r.items()))[1]) == set(r.values())

    i = 0
    for k, v in r:
        i += 1
        assert v in r.values()
        assert k in r.keys()
    assert i == 3

