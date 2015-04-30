from __future__ import unicode_literals, print_function, division
import numpy as np
import ahkab

class Test_DC_Solution:
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
        # setup analysis and simulate
        dca = ahkab.new_dc(-5, 5, 100, 'V1')
        self.r = ahkab.run(ttn, dca)['dc']

    def test(self):
        """Test results.dc_solution"""
        r = self.r
        # don't print this out, but still check we can convert
        # the results to string
        str(r)
        # is the independent variable correct?
        assert r.get_xlabel() == 'V1'
        assert r.units['V1'] == 'V'
        # and its values?
        assert np.alltrue(r.get_x() == np.array(list(ahkab.utilities.lin_axis_iterator(-5, 5, 100))))
        assert r['V1'][0] == -5
        assert np.allclose(r['V1'][-1], 5)
        assert (r.get_x() == r['V1']).all()
        # solution methods
        assert len(r) == len(r.variables)
        assert len(r) == len(list(r.keys()))
        assert r.has_key('V1')
        assert 'V1' in r
        assert not 'bogus' in r
        # Dictionary access with bogus keys
        try:
            r['sd']
            assert False
        except KeyError:
            pass
        # get() interface with default value
        assert r.get('sd', 'No such key') == 'No such key'
        assert not r.get('VN1').any() # all 0s
        # data shape
        assert r.asmatrix().shape == (8, 100)
        # values(), keys(), items()
        assert len(r.values()) == r.asmatrix().shape[0]
        assert len(r.values()[1]) == r.asmatrix().shape[1]
        assert set(list(zip(*r.items()))[0]) == set(r.keys())
        # Iterator interface
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



