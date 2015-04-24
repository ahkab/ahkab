from __future__ import print_function, division, unicode_literals

import numpy as np
import ahkab

from nose import SkipTest
from ahkab import py3compat

def test_pz_solution():
    """Test the PZ solution class"""
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

    print(str(r))

    np.allclose(r['p0']+r['p1'], -1034507*2, rtol=1e-3)
    np.allclose(abs(r['p0']-r['p1']), 107297253*2, rtol=1e-3)
    np.allclose(r['z0'], 0, rtol=1.)

    assert set(r.keys()) == {u'p0', u'p1', u'z0'}

    for i in {u'p0', u'p1', u'z0'}:
        assert i in r
        assert r.has_key(i)

    set(zip(*r.items())[0]) == {u'p0', u'p1', u'z0'}
    set(zip(*r.items())[1]) == set(r.values())

    i = 0
    for k, v in r:
        i += 1
        assert v in r.values()
        assert k in r.keys()
    assert i == 3

