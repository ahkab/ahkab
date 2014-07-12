from __future__ import print_function
import ahkab
import numpy as np
from ahkab.mosq import *

def mirror(vds, vgs, vbs):
    return -vds, vgs - vds, vbs - vds

def test_switch_nmos():
    mos_test_params = [dict(TYPE='n', KP=50e-6, VTO=.4)]
    mp = mos_test_params[0]
    m = mosq_mos_model(**mp)
    ma = mosq_device(1, 2, 3, 4, W=10e-6, L=1e-6, model=m)
    ma.part_id = "M1"

    for i in range(10):
        vds = np.random.rand()*8 - 4
        vgs = np.random.rand()*8 - 4
        vbs = np.random.rand()*8 - 4
        # print (vds, vgs, vbs)
        if vds < 0:
            assert (-vds, vgs - vds, vbs - vds) == \
                   m.get_voltages(vds, vgs, vbs)[0]
        else:
            assert (vds, vgs, vbs) == \
                   m.get_voltages(vds, vgs, vbs)[0]

def test_switch_pmos():
    mos_test_params = [dict(TYPE='p', KP=50e-6, VTO=-.4)]
    mp = mos_test_params[0]
    m = mosq_mos_model(**mp)
    ma = mosq_device(1, 2, 3, 4, W=10e-6, L=1e-6, model=m)
    ma.part_id = "M1"

    for i in range(10):
        vds = np.random.rand()*8 - 4
        vgs = np.random.rand()*8 - 4
        vbs = np.random.rand()*8 - 4
        print (vds, vgs, vbs)
        if vds > 0:
            assert (vds, - vgs + vds, - vbs + vds) == \
                   m.get_voltages(vds, vgs, vbs)[0]
        else:
            assert (-vds, -vgs, -vbs) == \
                   m.get_voltages(vds, vgs, vbs)[0]

def test_symm_nmos():
    mos_test_params = [dict(TYPE='n', KP=50e-6, VTO=.4)]
    mp = mos_test_params[0]
    m = mosq_mos_model(**mp)
    ma = mosq_device(1, 2, 3, 4, W=10e-6, L=1e-6, model=m)
    ma.part_id = "M1"
    for i in range(100):
        vds = np.random.rand()*8 - 4
        vgs = np.random.rand()*8 - 4
        vbs = np.random.rand()*8 - 4
        i1 = ma.istamp((vds, vgs, vbs))[1][0]
        vds2, vgs2, vbs2 = mirror(vds, vgs, vbs)
        i2 = ma.istamp((vds2, vgs2, vbs2))[1][0]
        assert i1 == -i2

def test_symm_pmos():
    mos_test_params = [dict(TYPE='p', KP=50e-6, VTO=-.4)]
    mp = mos_test_params[0]
    m = mosq_mos_model(**mp)
    ma = mosq_device(1, 2, 3, 4, W=10e-6, L=1e-6, model=m)
    ma.part_id = "M1"
    for i in range(100):
        vds = np.random.rand()*8 - 4
        vgs = np.random.rand()*8 - 4
        vbs = np.random.rand()*8 - 4
        i1 = ma.istamp((vds, vgs, vbs))[1][0]
        vds2, vgs2, vbs2 = mirror(vds, vgs, vbs)
        i2 = ma.istamp((vds2, vgs2, vbs2))[1][0]
        assert i1 == -i2

def test_symm_pn():
    mos_test_params = [dict(TYPE='p', KP=50e-6, VTO=-.4),
                       dict(TYPE='n', KP=50e-6, VTO=+.4)]
    mpp = mos_test_params[0]
    mpn = mos_test_params[1]
    mp = mosq_mos_model(**mpp)
    mn = mosq_mos_model(**mpn)
    mpa = mosq_device(1, 2, 3, 4, W=10e-6, L=1e-6, model=mp)
    mpa.part_id = "M1"
    mna = mosq_device(1, 2, 3, 4, W=10e-6, L=1e-6, model=mn)
    mna.part_id = "M2"
    for i in range(100):
        vds = np.random.rand()*8 - 4
        vgs = np.random.rand()*8 - 4
        vbs = np.random.rand()*8 - 4
        i1 = mna.istamp((vds, vgs, vbs))[1][0]
        vds2, vgs2, vbs2 = (-vds, -vgs, -vbs)
        i2 = mpa.istamp((vds2, vgs2, vbs2))[1][0]
        print(mn.get_voltages(vds, vgs, vbs)[0], '->', i1)
        print(mp.get_voltages(vds, vgs, vbs)[0], '->', i2)
        print(i1, i2)
        assert i1 == -i2

