from __future__ import unicode_literals, print_function, division
import numpy as np
import ahkab

class Test_TRAN_solution:
    def setUp(self):
        ttn = ahkab.Circuit('Twin-T Notch Stopband filter')
        rect_wave = ahkab.devices.pulse(v1=-1, v2=+1, td=0, tr=10e-6, pw=90e-6,
                                        tf=10e-6, per=200e-6)
        ttn.add_vsource('V1', 'in', ttn.gnd, dc_value=1, ac_value=1,
                        function=rect_wave)
        # first path
        ttn.add_capacitor('C1', 'in', 'n1', 2.2e-12)
        ttn.add_capacitor('C2', 'n1', 'out', 2.2e-12)
        ttn.add_resistor('R1', 'n1', ttn.gnd, 1e3)
        # second path
        ttn.add_resistor('R2', 'in', 'n2', 2e3)
        ttn.add_resistor('R3', 'n2', 'out', 2e3)
        ttn.add_capacitor('C3', 'n2', ttn.gnd, 2*2.2e-12)
        ttn.add_vcvs('E1', 'outb', ttn.gnd, 'out', ttn.gnd, 1.)
        # create a simulation object and run it!
        # using no step control, we know EXACTLY what the time
        # axis looks like. Easier for testing.
        tra = ahkab.new_tran(0, 50e-6, tstep=5e-6, x0=None,
                             use_step_control=False)
        self.r = ahkab.run(ttn, tra)['tran']

    def test(self):
        """Test results.tran_solution"""
        r = self.r
        # check no crashing
        str(r)
        # x-axis tests
        # the x0 is NOT added to the results.
        # hence no t=0
        assert np.alltrue(r.get_x() == np.array(list(ahkab.utilities.lin_axis_iterator(5e-6, 50e-6, 10))))
        assert r.get_xlabel() == 'T'
        assert (r.get_x() == r['T']).all()
        # solution methods
        assert r.asmatrix().shape == (8, 10)
        assert len(r) == len(r.variables)
        assert len(r) == len(list(r.keys()))
        assert r.has_key('T')
        assert 'T' in r
        assert not 'bogus' in r
        # x-axis values (we disabled the adaptive stepping)
        assert r['T'][0] == 5e-6
        assert r['T'][-1] >= 50e-6
        # dictionary access with bogus keys
        try:
            r['sd']
            assert False
        except KeyError:
            pass
        # default value in get()
        assert r.get('sd', 'No such key') == 'No such key'
        # keys, values, items
        assert len(r.values()) == r.asmatrix().shape[0]
        assert len(r.values()[1]) == r.asmatrix().shape[1]
        assert set(list(zip(*r.items()))[0]) == set(r.keys())
        # dictionary interface
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
                # no break occurred! something's off!
                assert False

