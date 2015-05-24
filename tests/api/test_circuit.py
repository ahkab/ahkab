# -*- coding: iso-8859-1 -*-
# test_circuit.py
# Utilities module test file for the Ahkab simulator
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
# along with ahkab.  If not, see <http://www.gnu.org/licenses/>.

"""
This module contains test functions for the circuit module.

"""

from __future__ import (unicode_literals, absolute_import,
                        division, print_function)

import ahkab
import copy

from nose.tools import raises

def test_remove_elem_linear():
    """Test circuit.remove_elem with linear elem"""
    c = ahkab.Circuit('test')
    c.add_resistor('R1', 'n1', 'n2', 1e3)
    c.add_resistor('R2', 'n2', c.gnd, 1e3)
    c.add_resistor('R3', 'n2', 'n3', 1e3)
    c.add_resistor('R4', 'n3', 'n4', 1e3)

    # save the internal node reference of n1 to check
    old_nodes_dict = copy.copy(c.nodes_dict)
    # remove R1, which should trigger the removal of n1
    c.remove_elem('R1')

    # check
    assert 'n1' not in c.nodes_dict
    assert old_nodes_dict['n1'] not in c.nodes_dict
    old_nodes_dict.pop(old_nodes_dict['n1'])
    old_nodes_dict.pop('n1')
    # nothing else got removed
    assert old_nodes_dict == c.nodes_dict
    assert 'R1' not in [e.part_id.upper() for e in c]

    # remove R3, which should trigger the removal of no
    # nodes
    # this time we call the method with the element
    R3 = [e for e in c if e.part_id.upper() == 'R3'][0]
    c.remove_elem(R3)
    # nothing got removed, right?
    assert old_nodes_dict == c.nodes_dict
    assert R3 not in c
    assert 'R3' not in [e.part_id for e in c]
    # YAY!


def test_find_vde_index():
    """Test circuit.find_vde_index()"""
    c = ahkab.Circuit('test')
    c.add_resistor('R1', 'n1', 'n2', 1e3)
    c.add_inductor('L1', 'n2', c.gnd, 1e-9)
    c.add_vsource('V1', 'n3', c.gnd, 1.)
    c.add_inductor('L2', 'n3', c.gnd, 1e-9)
    c.add_resistor('R2', 'n2', 'n3', 1e3)
    c.add_capacitor('C1', 'n3', 'n4', 1e-13)

    # test c.find_vde_index with elements as input
    Ls = [e for e in c if e.part_id[0].lower() == 'l']
    Ls.sort(key=lambda x: x.part_id)
    for L, i in zip(Ls, (0, 2)):
        assert c.find_vde_index(L) == i

    # test find_vde_index with IDs as input
    assert c.find_vde_index('l1') == 0
    assert c.find_vde_index('V1') == 1
    assert c.find_vde_index('L2') == 2


def test_find_vde():
    """Test circuit.find_vde()"""
    c = ahkab.Circuit('test')
    c.add_resistor('R1', 'n1', 'n2', 1e3)
    c.add_inductor('L1', 'n2', c.gnd, 1e-9)
    c.add_vsource('V1', 'n3', c.gnd, 1.)
    c.add_inductor('L2', 'n3', c.gnd, 1e-9)
    c.add_resistor('R2', 'n2', 'n3', 1e3)
    c.add_capacitor('C1', 'n3', 'n4', 1e-13)

    nnodes = len(c.nodes_dict)/2

    # test find_vde
    for pi in ('L1', 'L2', 'v1'):
        e = c.find_vde(nnodes - 1 + c.find_vde_index(pi))
        assert e.part_id.lower() == pi.lower()

@raises(ValueError)
def test_create_node1():
    """Test create_node() rejecting doubles"""
    cir = ahkab.Circuit('Circuit test')
    n1 = cir.create_node('n1')
    cir.add_capacitor('C1', 'n1', 'n2', 1e-6)
    cir.create_node('n1')

@raises(TypeError)
def test_create_node2():
    """Test create_node() rejecting nodes of wrong types"""
    cir = ahkab.Circuit('Circuit test')
    n1 = cir.create_node(1)

@raises(TypeError)
def test_add_node1():
    """Test add_node() rejecting nodes of wrong types"""
    cir = ahkab.Circuit('Circuit test')
    n1 = cir.add_node(1)

def test_add_node2():
    """Check the return value of add_node()"""
    cir = ahkab.Circuit('Circuit test')
    cir.add_node('n1')
    assert cir.nodes_dict['n1'] == cir.add_node('n1')

def test_add_node3():
    """Check add_node() correctly adding nodes 1/2"""
    cir = ahkab.Circuit('Circuit test')
    n1 = cir.add_node('n1')
    assert n1 in cir.nodes_dict

def test_add_node4():
    """Check add_node() correctly adding nodes 2/2"""
    cir = ahkab.Circuit('Circuit test')
    n1 = cir.add_node('n1')
    assert cir.nodes_dict[n1] in cir.nodes_dict

def test_create_node3():
    """Check create_node() correctly adding nodes 1/2"""
    cir = ahkab.Circuit('Circuit test')
    n1 = cir.create_node('n1')
    assert n1 in cir.nodes_dict

def test_add_node4():
    """Check create_node() correctly adding nodes 2/2"""
    cir = ahkab.Circuit('Circuit test')
    n1 = cir.create_node('n1')
    assert cir.nodes_dict[n1] in cir.nodes_dict

def test_create_node5():
    """Check the return value of create_node()"""
    cir = ahkab.Circuit('Circuit test')
    assert 'n1' == cir.create_node('n1')

@raises(TypeError)
def test_new_internal_node1():
    """Check the internal nodes exceptions"""
    cir = ahkab.Circuit('Circuit test')
    n = cir.new_internal_node()
    assert cir.is_int_node_internal_only(n)

def test_new_internal_node2():
    """Check the internal nodes"""
    cir = ahkab.Circuit('Circuit test')
    n = cir.new_internal_node()
    assert cir.is_int_node_internal_only(cir.nodes_dict[n])

def test_get_nodes_number():
    """Test get_nodes_number()"""
    cir = ahkab.Circuit('Circuit test')
    n1 = cir.create_node('n1')
    n2 = cir.create_node('n2')
    n3 = cir.create_node('n3')
    n4 = cir.create_node('n4')
    n5 = cir.create_node('n5')
    n1 = cir.add_node('n1')
    gnd = cir.create_node('0')
    assert cir.get_nodes_number() == 6
    #assert cir.nodes_dict[n1] in cir.nodes_dict
    assert cir.get_nodes_number() == int(len(cir.nodes_dict)/2)

if __name__ == '__main__':
    test_remove_elem_linear()
    test_find_vde()
    test_find_vde_index()
    test_create_node1()
    test_create_node2()
    test_add_node1()
    test_add_node2()
    test_add_node3()
    test_add_node4()
    test_create_node3()
    test_add_node4()
    test_create_node5()
    test_new_internal_node1()
    test_new_internal_node2()
    test_get_nodes_number()
