import copy
import sys
import os
import numpy as np

from nose.plugins.skip import SkipTest

from ahkab import main

def get_reference_path():
    test_id = 'pz_spice'
    wd = os.getcwd()
    if os.path.split(wd)[1] == 'ahkab':
        reference_path = os.path.join(wd, 'tests', test_id)
    elif os.path.split(wd)[1] == 'tests':
        reference_path = os.path.join(wd, test_id)
    else:
        reference_path = "."
    return reference_path

def _check_singularities(res, ref, atol=1e-4, rtol=1e-3):
    ref = copy.copy(ref)
    for i in res:
        for si, s in enumerate(ref):
            if np.allclose(i, s, atol=atol, rtol=rtol):
                break
            assert si != len(ref) - 1
        ref.pop(si)
    assert not len(ref)

def test_pz_spice1():
    """PZ Test: 5TH-ORDER LOW_PASS FILTER"""
    netlist_path = os.path.join(get_reference_path(), 'pz_spice1.ckt')
    r = main(netlist_path, verbose=0, outfile='tmp')
    res = r['pz']
    poles = [res[s] for s in res.keys() if s[0] == 'p']
    poles_ref = [-1.105884e-02 -7.435365e-02j, 
                 -1.105884e-02 +7.435365e-02j,
                 -1.882392e-02 -1.418852e-01j,
                 -1.882392e-02 +1.418852e-01j,
                 -8.675361e-02+0.000000e+00j]
    _check_singularities(poles, poles_ref, atol=1e-4, rtol=1e-3)

    zeros = [res[s] for s in res.keys() if s[0] == 'z']
    zeros_ref = [-2.047019e-01j,
                 -8.164476e-02j,
                 8.164476e-02j,
                 2.047019e-01j]
    _check_singularities(zeros, zeros_ref, atol=1e-4, rtol=1e-3)

def test_pz_spice2():
    """PZ Test: KERWIN'S CIRCUIT - JW-AXIS TRANSMISSION ZEROS"""
    netlist_path = os.path.join(get_reference_path(), 'pz_spice2.ckt')
    r = main(netlist_path, verbose=0, outfile='tmp')
    res = r['pz']
    poles = [res[s] for s in res.keys() if s[0] == 'p']
    poles_ref = [-7.964016e-03 + 1.589515e-01j,
                 -7.964016e-03 -1.589515e-01j,
                 -2.250812e-01 + 0.000000e+00j]
    _check_singularities(poles, poles_ref, atol=1e-4, rtol=1e-3)

    zeros = [res[s] for s in res.keys() if s[0] == 'z']
    zeros_ref = [-2.250812e-01j, +2.250812e-01j, -2.250812e-01]
    _check_singularities(zeros, zeros_ref, atol=1e-4, rtol=1e-3)

def test_pz_spice3():
    """PZ Test: High-Pass Butterworth Filter"""
    if 'TRAVIS' in os.environ:
        raise SkipTest
    netlist_path = os.path.join(get_reference_path(), 'pz_spice3.ckt')
    r = main(netlist_path, verbose=0, outfile='tmp')
    res = r['pz']
    poles = [res[s] for s in res.keys() if s[0] == 'p']
    poles_ref = [-6.090889e-02 + 1.470617e-01j,
                 -6.090889e-02 - 1.470617e-01j,
                 -1.470254e-01 + 6.093849e-02j,
                 -1.470254e-01 - 6.093849e-02j]

    _check_singularities(poles, poles_ref, atol=1e-4, rtol=1e-3)

    zeros = [res[s] for s in res.keys() if s[0] == 'z']
    zeros_ref = [0+0j, 0+0j, 0+0j, 0+0j]
    _check_singularities(zeros, zeros_ref, atol=1e-4, rtol=1e-3)

def test_pz_spice4():
    """PZ Test: Simple Amplifier"""
    netlist_path = os.path.join(get_reference_path(), 'pz_spice4.ckt')
    r = main(netlist_path, verbose=0, outfile='tmp')
    res = r['pz']
    poles = [res[s] for s in res.keys() if s[0] == 'p']
    poles_ref = [-2.248151e+05 +0j, -2.253434e+07 +0j]
    _check_singularities(poles, poles_ref, atol=1e2, rtol=1e-3)

    zeros = [res[s] for s in res.keys() if s[0] == 'z']
    zeros_ref = [6.366198e+08+0j]
    _check_singularities(zeros, zeros_ref, atol=1e2, rtol=1e-3)

