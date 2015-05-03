from __future__ import unicode_literals, print_function, division
import numpy as np
import ahkab

class Test_PSS_solution:
    def setUp(self):
        ttn = ahkab.Circuit('PSS Linear circuit')
        sin_wave = ahkab.devices.sin(vo=0, va=2, td=0, freq=1e6, theta=0, phi=0.)
        ttn.add_vsource('VIN', 'in', ttn.gnd, dc_value=5, function=sin_wave)
        ttn.add_resistor('R1', 'in', 'n1', 10e3)
        ttn.add_resistor('R2', 'n1', 'out', 20e3)
        ttn.add_resistor('R3', 'n2', ttn.gnd, 10e3)
        ttn.add_resistor('R4', 'n3', ttn.gnd, 20e3)
        ttn.add_resistor('R5', 'n3', 'out', 40e3)
        ttn.add_capacitor('C1', 'n1', 'n2', 31.83e-12)
        ttn.add_capacitor('C2', 'n1', 'n2', 15.91e-12)
        ttn.add_vcvs('E1', 'out', ttn.gnd, 'n2', 'n3', 1e6)
                # create a simulation object and run it!
        op = ahkab.new_op()
        pssa = ahkab.new_pss(1e-6, points=60)
        self.r = ahkab.run(ttn, [op, pssa])['pss']

    def test(self):
        """Test results.pss_solution"""
        r = self.r
        # check no crashing
        str(r)
        # x-axis tests
        assert np.allclose(r.get_x(), np.array(list(ahkab.utilities.lin_axis_iterator(0, 1e-6, 60))))
        assert r.get_xlabel() == 'T'
        # solution methods
        assert len(r) == len(list(r.keys()))
        assert r.has_key('VIN')
        assert 'VIN' in r
        assert not 'bogus' in r
        # units check
        assert r.units['T'] == 's'
        assert r.units['ViN'] == 'V'
        # asarray()
        assert r.asarray().shape == (8, 60)
        assert np.alltrue(r.asarray()[0, :] == r.get_x())
        # dictionary access with bogus cells
        try:
            r['sd']
            assert False
        except KeyError:
            pass
        # get()
        assert r.get('sd', 'No such key') == 'No such key'
        # keys, values, items
        assert len(list(r.keys())) == len(list(r.values()))
        assert set(list(zip(*r.items()))[0]) == set(list(r.keys()))
        # iterator
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

