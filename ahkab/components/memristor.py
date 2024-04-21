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
def g(self):
    return self._value


