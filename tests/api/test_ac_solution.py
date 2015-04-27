from __future__ import unicode_literals, print_function, division
import numpy as np
import ahkab


class Test_AC_solution:
    def setUp(self):
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
        aca = ahkab.new_ac(1e7, 1e10, 100, x0=None)
        self.r = ahkab.run(ttn, aca)['ac']

    def test(self):
        """Test results.ac_solution"""
        r = self.r
        str(r)
        # keyerrors
        try:
            r['sd']
            assert False
        except KeyError:
            pass
        #get method
        assert r.get('sd', 'No such key') == 'No such key'
        assert r.get('VN1').any()
        # axis endpoints
        assert r['w'][0] == 10000000
        assert np.allclose(r['w'][-1], 10000000000)
        # data shape
        assert r.asmatrix().shape == (8, 100)
        # data
        assert len(r.values()) == r.asmatrix().shape[0]
        assert len(r.values()[1]) == r.asmatrix().shape[1]
        # iterator checks
        keys = []
        values = []
        for k, v in r:
            keys.append(k)
            values.append(v)
        assert len(keys) == len(r.keys())
        assert len(values) == len(r.values())
        for v in r.values():
            for vi in values:
                if np.allclose(vi, v):
                    break
            else:
                assert False
        # items, values, keys etc
        assert set(list(zip(*r.items()))[0]) == set(r.keys())

