# -*- coding: iso-8859-1 -*-
# Copyright 2006 Giuseppe Venturini

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
# along with ahkab.  If not, see <http://www.gnu.org/licenses/>.
from .Component import Component
import numpy as np
class Memristor(Component):
    """
    A memristor

    """
    
    def __init__(self, part_id, n1, n2, value):
        self.part_id = part_id
        self._value = value
        self.is_nonlinear = True
        self.is_symbolic = True
        self.n1 = n1
        self.n2 = n2

    @property
    def g(self,v=0,value = 0):
        return self._value

    @g.setter
    def g(self,g):
        self._value = g
    
    def get_ports(self):
        return (self.n1,self.n2)

    def i(self,v,time):
        return self._value*v

    def _sig(self, v):
        return 1/(1 + np.exp(-v))

    def EulerForward(self, dt, v):
        tau = 1
        dgdt = (self._sig(v) - self._value)/tau
        g_new =+ dt*dgdt
        return g_new

    def get_netlist_elem_line(self, nodes_dict):
        return "%s %s %s %g" % (self.part_id, nodes_dict[self.n1],
                                nodes_dict[self.n2], self.value)


