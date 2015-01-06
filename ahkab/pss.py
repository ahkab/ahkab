# -*- coding: iso-8859-1 -*-
# pss.py
# Periodic Steady State analysis
# Copyright 2011-2013 Giuseppe Venturini

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

from __future__ import (unicode_literals, absolute_import,
                        division, print_function)

from . import shooting
from . import bfpss
from . import options

# Periodic Steady State (PSS) analysis module.

# This module is an interface to a generic PSS analysis,
# which will be set up automatically for you according to
# the algorithm selected.

# Netlist syntax
# .PSS PERIOD=n [points=n step=n autonomous=bool method=str]

specs = {'pss': {'tokens': ({
                            'label': 'period',
                            'pos': 0,
                            'type': float,
                            'needed': True,
                            'dest': 'period',
                            'default': None
                            },
        {
                            'label': 'points',
                            'pos': None,
                            'type': float,
                            'needed': True,
                            'dest': 'points',
                            'default': None
                            },
                 {
                            'label': 'step',
                            'pos': None,
                            'type': float,
                            'needed': False,
                            'dest': 'step',
                            'default': None
                            },
                 {
                            'label': 'autonomous',
                            'pos': None,
                            'type': bool,
                            'needed': False,
                            'dest': 'autonomous',
                            'default': '0'
                            },
                 {
                            'label': 'method',
                            'pos': None,
                            'type': str,
                            'needed': False,
                            'dest': 'method',
                            'default': options.SHOOTINGPSS
                            }
                            )
                 }
         }


def pss_analysis(*largs, **args):
    """Perform a PSS analysis.

    The only required argument is ``method`` (string), which
    selects the algorithm to be used.

    Two algorithms are a shooting PSS, selected with the ``"shooting"`` switch
    and brute-force PSS, selected by the ``"brute-force"`` switch. Any other
    value for the ``method`` parameter will result in a ``ValueError`` exception
    being raised.

    The rest of the arguments will be passed to the algorithm implementing the
    analysis.

    For the shooting algorithm see :mod:`ahkab.shooting`, for the brute-force
    algorithm see :mod:`ahkab.bfpss`.

    **Returns:**
    
    sol : PSS solution object (:class:`results.pss_solution`)
        The solution.
    """
    m = args.pop('method').lower()
    if m == options.SHOOTINGPSS:
        r = shooting.shooting_analysis(*largs, **args)
    elif m == options.BFPSS:
        r = bfpss.bfpss(*largs, **args)
    else:
        raise Exception("Unknown PSS method %s" % m)
    return r
