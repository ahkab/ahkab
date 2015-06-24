#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-
# ahkab.py
# The frontend of the simulator
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

""" ahkab is an easy electronic circuit simulator.
"""

from __future__ import print_function

import os

# Package testing can be done remotely, without display. This would make
# matplotlib fail (and, sometimes, the test as well).
# If we check for $DISPLAY, that makes us probably lose in portability,
# because does Windows have the DISPLAY env variable defined?
try:
    import matplotlib
    if os.system('python -c "import matplotlib.pyplot as plt;plt.figure()"'):
        # no display, use a display-less backend
        matplotlib.use('Agg')
except ImportError:
    # plotting.py will complain about this later on.
    pass

from .ahkab import new_op, new_dc, new_tran, new_ac, new_pss, new_pz
from .ahkab import new_symbolic, queue, run, new_x0, icmodified_x0
from .ahkab import get_op_x0, set_temperature, process_postproc, main
from .__version__ import __version__
from .circuit import Circuit

__all__ = ['new_op', 'new_dc', 'new_tran', 'new_ac', 'new_pss',
           'new_symbolic', 'new_pz', 'queue', 'run', 'new_x0',
           'get_op_x0', 'set_temperature', 'main', 'Circuit']
