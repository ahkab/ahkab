#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-
# __main__.py
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

from __future__ import print_function, division

import sys
from optparse import OptionParser

from . import ahkab
from . import options
from . import transient
from . import utilities
from .__version__ import __version__

def _cli():
    usage = "usage: \t%prog [options] <filename>\n\nThe filename is the " + \
            "netlist to be open. Use - (a dash) to read from stdin."
    parser = OptionParser(usage, version="%prog " + __version__ + \
             " (c) 2006-2015 Giuseppe Venturini")

    # general options
    parser.add_option("-v", "--verbose", action="store", type="string",
                     dest="verbose", default="3", help="Verbose level: " +
                    "from 0 (almost silent) to 5 (debug)")
    parser.add_option("-p", "--print", action="store_true",
                      dest="print_circuit", default=False, help="Print " +
                      "the parsed circuit")
    parser.add_option("-o", "--outfile", action="store", type="string",
                      dest="outfile", default='stdout', help="Data output " +
                      "file. Defaults to stdout.")
    parser.add_option("", "--dc-guess", action="store", type="string",
                      dest="dc_guess", default="guess", help="Guess to be " +
                      "used to start a OP or DC analysis: none or guess. " +
                      "Defaults to guess.")
    parser.add_option("-t", "--tran-method", action="store", type="string",
                      dest="method", default=transient.TRAP.lower(),
                      help="Method to be used in transient analysis: " +
                      transient.IMPLICIT_EULER.lower() + ", " +
                      transient.TRAP.lower() + ", " +
                      transient.GEAR2.lower() + ", " +
                      transient.GEAR3.lower() + ", " +
                      transient.GEAR4.lower() + ", " +
                      transient.GEAR5.lower() + " or " +
                      transient.GEAR6.lower() + ". Defaults to TRAP.")
    parser.add_option("", "--t-fixed-step", action="store_true",
                      dest="no_step_control", default=False, help="Disables" +
                      " the step control in transient analysis.")
    parser.add_option("", "--v-absolute-tolerance", action="store",
                      type="string", dest="vea", default=None, help="Voltage" +
                      " absolute tolerance. Default: " + str(options.vea) + " V")
    parser.add_option("", "--v-relative-tolerance", action="store",
                      type="string", dest="ver", default=None, help="Voltage " +
                      "relative tolerance. Default: " + str(options.ver))
    parser.add_option("", "--i-absolute-tolerance", action="store",
                      type="string", dest="iea", default=None, help="Current " +
                      "absolute tolerance. Default: " + str(options.iea) + " A")
    parser.add_option("", "--i-relative-tolerance", action="store",
                      type="string", dest="ier", default=None, help="Current " +
                      "relative tolerance. Default: " + str(options.ier))
    parser.add_option("", "--h-min", action="store", type="string",
                      dest="hmin", default=None, help="Minimum time step. " +
                      "Default: " + str(options.hmin))
    parser.add_option("", "--dc-max-nr", action="store", type="string",
                      dest="dc_max_nr_iter", default=None, help="Maximum " +
                      "number of NR iterations for DC and OP analyses. " +
                      "Default: " + str(options.dc_max_nr_iter))
    parser.add_option("", "--t-max-nr", action="store", type="string",
                      dest="transient_max_nr_iter", default=None,
                      help="Maximum number of NR iterations for each time " +
                      "step during transient analysis. Default: " +
                      str(options.transient_max_nr_iter))
    parser.add_option("", "--t-max-time", action="store", type="string",
                      dest="transient_max_time_iter", default=None,
                      help="Maximum number of time iterations during " +
                      "transient analysis. Setting it to 0 (zero) disables " +
                      "the limit. Default: " +
                      str(options.transient_max_time_iter))
    parser.add_option("", "--s-max-nr", action="store", type="string",
                      dest="shooting_max_nr_iter", default=None,
                      help="Maximum number of NR iterations during shooting " +
                      "analysis. Setting it to 0 (zero) disables the " +
                      "limit. Default: " + str(options.shooting_max_nr_iter))
    parser.add_option("", "--gmin", action="store", type="string", dest="gmin",
                      default=None, help="The minimum conductance to ground. " +
                      "Inserted when requested. Default: " +
                      str(options.gmin))
    parser.add_option("", "--cmin", action="store", type="string", dest="cmin",
                      default=None, help="The minimum capacitance to ground. " +
                      "Default: " + str(options.cmin))

    (cli_options, remaning_args) = parser.parse_args()

    verbose = int(cli_options.verbose)
    if cli_options.method is not None:
        options.default_tran_method = cli_options.method.upper()
    if cli_options.vea is not None:
        options.vea = float(cli_options.vea)
    if cli_options.ver is not None:
        options.ver = float(cli_options.ver)
    if cli_options.iea is not None:
        options.iea = float(cli_options.iea)
    if cli_options.ier is not None:
        options.ier = float(cli_options.ier)
    if cli_options.hmin is not None:
        options.hmin = float(cli_options.hmin)
    if cli_options.dc_max_nr_iter is not None:
        options.dc_max_nr_iter = int(float(cli_options.dc_max_nr_iter))
    if cli_options.transient_max_nr_iter is not None:
        options.transient_max_nr_iter = int(
            float(cli_options.transient_max_nr_iter))
    if cli_options.transient_max_time_iter is not None:
        options.transient_max_time_iter = int(
            float(cli_options.transient_max_time_iter))
    if cli_options.shooting_max_nr_iter is not None:
        options.shooting_max_nr_iter = int(
            float(cli_options.shooting_max_nr_iter))
    if cli_options.gmin is not None:
        options.gmin = float(cli_options.gmin)
    if cli_options.cmin is not None:
        options.cmin = float(cli_options.cmin)

    if not len(remaning_args) == 1:
        print("Usage: ./ahkab.py [options] <filename>\n./ahkab.py -h for help")
        sys.exit(1)
    if remaning_args[0] == '-':
        read_netlist_from_stdin = True
    else:
        read_netlist_from_stdin = False
    if not read_netlist_from_stdin and not utilities.check_file(remaning_args[0]):
        sys.exit(23)

    options.transient_no_step_control = cli_options.no_step_control
    options.dc_use_guess = cli_options.dc_guess
    options.cli = True
    ahkab._print = cli_options.print_circuit

    # Program execution
    ahkab.main(filename=remaning_args[0],
               outfile=cli_options.outfile, verbose=verbose)

    sys.exit(0)

if __name__ == "__main__":
    _cli()
