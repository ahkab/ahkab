# -*- coding: iso-8859-1 -*-
# test_cid.py
# Unit tests for the case insensitive dictionary class
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
import numpy as np
import ahkab


class TestCID:
    def setUp(self):
        self.a = ahkab.results.case_insensitive_dict()
        self.a.update({'a':1, 'b':2})
        self.a.update({'c':3})


    def test(self):
        """Test results.case_insensitive_dict"""
        # str
        print(str(self.a))
        assert len(self.a) == 3
        # we said case insensitive!
        assert 'a' in self.a
        assert 'B' in self.a
        # basic stuff
        assert len(self.a) == len(list(self.a.keys()))
        assert self.a.has_key('C')
        assert not 'bogus' in self.a
        # key errors
        try:
            self.a['sd']
            assert False
        except KeyError:
            pass
        # get
        assert self.a.get('a') == self.a['a']
        # default values in get
        assert self.a.get('sd', 'No such key') == 'No such key'
        # keys, values, items
        assert len(list(self.a.keys())) == len(list(self.a.values()))
        assert set(list(zip(*self.a.items()))[0]) == set(list(self.a.keys()))
        assert set(list(zip(*self.a.items()))[1]) == set(list(self.a.values()))
        # iterator
        keys = []
        values = []
        for k in self.a:
            keys.append(k)
            values.append(self.a[k])
        assert len(keys) == len(self.a.keys())
        assert len(values) == len(self.a.values())
        for v in self.a.values():
            for vi in values:
                if vi == v:
                    break
            else:
                assert False

