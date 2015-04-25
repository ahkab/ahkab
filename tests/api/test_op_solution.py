from __future__ import unicode_literals, print_function, division

import numpy as np

import ahkab

def test_op_solution():
    """Test results.op_solution"""
    ######### CIRCUIT ##############
    ttn = ahkab.Circuit('Twin-T Notch Stopband filter')
    ttn.add_vsource('V1', 'in', ttn.gnd, dc_value=1, ac_value=1)
    # first path
    ttn.add_capacitor('C1', 'in', 'n1', 2.2e-12)
    ttn.add_capacitor('C2', 'n1', 'out', 2.2e-12)
    ttn.add_resistor('R1', 'n1', ttn.gnd, 1e3)
    # second path
    ttn.add_resistor('R2', 'in', 'n2', 2e3)
    ttn.add_resistor('R3', 'n2', 'out', 2e3)
    ttn.add_capacitor('C3', 'n2', ttn.gnd, 2*2.2e-12)
    ttn.add_vcvs('E1', 'outb', ttn.gnd, 'out', ttn.gnd, 1.)
    # set up the OP and run
    opa = ahkab.new_op()
    r = ahkab.run(ttn, opa)['op']
    
    ####### CHECKS ###########
    # str representation
    print(str(r))
    r.keys()

    # 'sd' is not an existing key
    try:
        r['sd']
        assert False
    except KeyError:
        pass

    # fallback on default
    assert r.get('sd', 1e3) == 1e3
    assert r.get('VN1') == 0

    # the important part is not the value
    np.allclose(r.asmatrix(), 
                np.array([[  1.00000000e+00],
                         [  0.00000000e+00],
                         [  1.00000000e+00],
                         [  1.00000000e+00],
                         [  1.00000000e+00],
                         [  1.01736171e-20],
                         [  0.00000000e+00]]), rtol=1e-3)

    r.print_short()
    set(list(zip(*r.items()))[0]) == set(r.keys())
    set(list(zip(*r.items()))[1]) == set(r.values())

    # iterator
    keys = set()
    values = set()
    for k, v in r:
        keys |= {k}
        values |= {float(v)}
    assert keys == set(r.keys())
    assert values == set(r.values())

