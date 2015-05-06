# -*- coding: iso-8859-1 -*-
# test_symbolic_solution.py
# Unit tests for the symbolic solution class
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

from __future__ import unicode_literals, print_function, division
import os
import numpy as np
import ahkab

class Test_symbolic_solution(object):
    def setUp(self):
        ttn = ahkab.Circuit('Twin-T Notch Stopband filter')
        ttn.add_vsource('V1', 'in', ttn.gnd, dc_value=1)
        # first path
        ttn.add_capacitor('C1', 'in', 'n1', 2.2e-12)
        ttn.add_capacitor('C2', 'n1', 'out', 2.2e-12)
        ttn.add_resistor('R1', 'n1', ttn.gnd, 1e3)
        # second path
        ttn.add_resistor('R2', 'in', 'n2', 2e3)
        ttn.add_resistor('R3', 'n2', 'out', 2e3)
        ttn.add_capacitor('C3', 'n2', ttn.gnd, 2*2.2e-12)
        ttn.add_vcvs('E1', 'outb', ttn.gnd, 'out', ttn.gnd, 1.)
        sa = ahkab.new_symbolic(source=None)
        self.r = ahkab.run(ttn, sa)['symbolic'][0]

    def test(self):
        """Test results.symbolic_solution"""
        r = self.r
        # only check no choking on str()
        str(r)
        # dictionary access
        assert len(r) == len(list(r.keys()))
        assert r.has_key('VIN')
        assert 'VIN' in r
        assert not 'bogus' in r
        # symbol-related functionality
        assert r.as_symbol('R1') == r.as_symbols('R1')[0]
        assert str(r.as_symbol('R1')) == 'R1'
        assert [r.as_symbol('R1'), r.as_symbol('R2'), r.as_symbol('R3')] == \
                r.as_symbols('R1 R2 R3')
        # load and save
        r.filename = './tmp_test_file_for_symbolic_results.pickle'
        r.save()
        s = ahkab.results.symbolic_solution.load(r.filename)
        # dictionary interface with bogus keys
        try:
            r['sd']
            assert False
        except KeyError:
            pass
        # get() and its default value
        assert r.get('sd', 'No such key') == 'No such key'
        # keys, values and items
        assert len(list(r.keys())) == len(list(r.values()))
        assert set(list(zip(*r.items()))[0]) == set(list(r.keys()))
        assert set(list(zip(*r.items()))[1]) == set(list(r.values()))
        # iterator test
        keys = []
        values = []
        for k, v in r:
            keys.append(k)
            values.append(v)
        assert len(keys) == len(list(r.keys()))
        assert len(values) == len(list(r.values()))
        for v in r.values():
            for vi in values:
                if vi == v:
                    break
            else:
                assert False
        os.remove(r.filename)

