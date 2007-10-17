# -*- coding: iso-8859-1 -*-
# options.py
# Global options file
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

# global: errors
# voltage: absolute and relative tolerances
vea = 1e-6
ver = 1e-3
# current: absolute and relative tolerances
iea = 1e-9
ier = 1e-3

# global: circuit
gmin = 1e-12
print_int_nodes = True

# global: solving
nr_dump_first_iters = False  	# the first iterations will be dumped see dc_analysis.get_td()
nl_voltages_lock = True 	# Apply damping - slows down solution.
nl_voltages_lock_factor = 4  	# if we allow the voltage on controlling ports to change too much,
				# we may have current overflows. Think about a diode (exp).
				# So we allow them to change of nl_voltages_lock_factor*Vth at most
				# and damp all variables accordingly

# dc
dc_max_nr_iter = 10000
dc_max_guess_effort = 250000

# transient
default_tran_method = "IMPLICIT_EULER"
hmin = 1e-20
transient_max_time_iter = 0 #disabled
transient_max_nr_iter = 20
# use the prediction value as first guess for x(n+1), otherwise uses x(n)
transient_prediction_as_x0 = True
transient_use_aposteriori_step_control = True
# we do not want to redo the iteraction if the aposteriori check suggests a step that is 
# very close to the one we already used. 0.9 seems to be a good idea.
transient_aposteriori_step_threshold = 0.9

# shooting
shooting_default_points = 100
shooting_max_nr_iter = 10000
