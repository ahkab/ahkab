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

class Mysistor(Component):
    """A memristor.

    .. image:: images/elem/resistor.svg

    **Parameters:**

    part_id : string
        The unique identifier of this element. The first letter should be
        ``'R'``.
    n1 : int
        *Internal* node to be connected to the anode.
    n2 : int
        *Internal* node to be connected to the cathode.
    value : float
        Resistance in ohms.

     """
    #
    #             /\    /\    /\
    #     n1 o---+  \  /  \  /  \  +---o n2
    #                \/    \/    \/
    #
    def __init__(self, part_id, n1, n2, value, rho_b=0.1, length_channel=10e-6, radius_base=200e-9, tau=0.0048, pressure=0, delta_rho=0, has_converged=False):


        # ID and nodes
        self.part_id = part_id
        self.is_nonlinear = False
        self.is_symbolic = True
        self.n1 = n1
        self.n2 = n2

        # convergence
        self.has_converged = has_converged
        
        # read inputs
        self._value = value     # reiststance [Ohm]
        self._g = 1./value      # conductance
        self.tau = tau      # relaxation period [s]   
        N_A = 6.022e23      #Avogadro numeber   [1]
        # print(rho_b)
        self.rho_b = rho_b * N_A     #bulk density [m^-3]
        self.length_channel = length_channel 
        self.initial_length = length_channel
        self.radius_base = radius_base
        self.pressure = pressure
        self.delta_rho = delta_rho
        # print(self.delta_rho)
        
        self.radius_tip = 50e-9     # [m]
        # self.radius_base = self.rbrt*self.radius_tip
        # self.radius_base = self.rb
        self.delta_radius = self.radius_base - self.radius_tip
        average_radius = (self.radius_base**2 + self.radius_tip**2 + self.radius_base*self.radius_tip)/3


        # constants
        electron_charge = 1.602e-19   # [C] = [A s]
        kbT = (1.38e-23)*(293.15)      # [J] = [kg m^2 s^-2]
        eta = 1.01e-3      #viscosity [mPa s] = [kg m^-1 s^-1] 
        epsilon = 0.71e-9     # [F m^-1] = [kg^-1 m^-1 s^4 A^2]
        sigma = -0.0015e18     #surface charge [m^-2]
        phi0 = -10e-3      # [V] = [kg m^2 s^-3 A^-1]
        diff_coefficient = 1.75e-9      # [m^2 s^-1]

        # derived constants
        w = (electron_charge*diff_coefficient*eta)/(kbT*epsilon*phi0)
        Du = sigma/(2*self.rho_b*self.radius_tip)

        # Ohmic conductance
        g_1 = np.pi*self.radius_tip*self.radius_base/self.length_channel
        g_2 = 2*self.rho_b*(electron_charge**2)*diff_coefficient/kbT
        
        self.g_0 = g_1*g_2*1e12     # [S]

        # ONLY VOLTAGE
        self.delta_g = (-2)*(w)*self.delta_radius*Du/self.radius_base
        
        # VOLTAGE AND PRESSURE
        
        self.density_inhomo_over_potential = (2*(self.delta_radius)*sigma*electron_charge)/(kbT*self.radius_tip**2)
        
        self.q_potential = (-np.pi)*self.radius_base*self.radius_tip*epsilon*phi0/(eta*self.length_channel)
        
        self.q_pressure = (np.pi*(self.radius_base*self.radius_tip)**3)/(8*self.length_channel*eta*average_radius)

        self.peclet_over_q = length_channel/(diff_coefficient*np.pi*self.radius_tip**2)


    def init_derived_param(self):

        # constants
        electron_charge = 1.602e-19   # [C] = [A s]
        kbT = (1.38e-23)*(293.15)      # [J] = [kg m^2 s^-2]
        eta = 1.01e-3      #viscosity [mPa s] = [kg m^-1 s^-1] 
        epsilon = 0.71e-9     # [F m^-1] = [kg^-1 m^-1 s^4 A^2]
        sigma = -0.0015e18     #surface charge [m^-2]
        phi0 = -10e-3      # [V] = [kg m^2 s^-3 A^-1]
        diff_coefficient = 1.75e-9      # [m^2 s^-1]

        # derived constants
        w = (electron_charge*diff_coefficient*eta)/(kbT*epsilon*phi0)
        Du = sigma/(2*self.rho_b*self.radius_tip)
        average_radius = (self.radius_base**2 + self.radius_tip**2 + self.radius_base*self.radius_tip)/3

        # Ohmic conductance
        g_1 = np.pi*self.radius_tip*self.radius_base/self.length_channel
        g_2 = 2*self.rho_b*(electron_charge**2)*diff_coefficient/kbT
        
        self.g_0 = g_1*g_2*1e12     # [S]

        # ONLY VOLTAGE
        self.delta_g = (-2)*(w)*self.delta_radius*Du/self.radius_base
        
        # VOLTAGE AND PRESSURE
        
        self.density_inhomo_over_potential = (2*(self.delta_radius)*sigma*electron_charge)/(kbT*self.radius_tip**2)
        
        self.q_potential = (-np.pi)*self.radius_base*self.radius_tip*epsilon*phi0/(eta*self.length_channel)
        
        self.q_pressure = (np.pi*(self.radius_base*self.radius_tip)**3)/(8*self.length_channel*eta*average_radius)

        self.peclet_over_q = self.length_channel/(diff_coefficient*np.pi*self.radius_tip**2)
        

    @property
    def g(self, v=0, time=0):
        return self._g

    @g.setter
    def g(self, g):
        self._g = g
        self._value = 1./g

    @property
    def value(self, v=0, time=0):
        return self._value

    @value.setter
    def value(self, value):
        self._value = value
        self._g = 1./value

    def new_nodes(self, n1, n2):
        self.n1 = n1
        self.n2 = n2

    def i(self, v, time=0):
        return 0

    def get_op_info(self, ports_v):
        """Information regarding the Operating Point (OP)

        **Parameters:**

        ports_v : list of lists
            The parameter is to be set to ``[[v]]``, where ``v`` is the voltage
            applied to the resistor terminals.

        **Returns:**

        op_keys : list of strings
            The labels corresponding to the numeric values in ``op_info``.
        op_info : list of floats
            The values corresponding to ``op_keys``.
        """
        vn1n2 = float(ports_v[0][0].item())
        in1n2 = float(ports_v[0][0].item()/self.value)
        power = float(ports_v[0][0].item() ** 2 / self.value)
        
        op_keys = ['Part ID', u"R [\u2126]", "V(n1,n2) [V]", "I(n1->n2) [A]", "P [W]"]
        op_info = [self.part_id.upper(), self.value, vn1n2, in1n2, power]
        return op_keys, op_info
