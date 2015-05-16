Command line help
=================

The ``ahkab`` simulator has a command line interface that
allows for quick simulation of netlist decks, without the need
to load the Python interpreter explicitely.

Several switches are available, to set the input and output files
and to override some built-in options.

Notice that options set on the command line always take precedence
on any netlist option or any value set in :mod:`ahkab.options`.

Usage:
------

::

     ahkab [options] <filename>

The filename is the netlist to be open. Use - (a dash) to read from stdin.

Options:
--------

--version 
    show program's version number and exit
-h, --help 
    show this help message and exit
-v VERBOSE, --verbose=VERBOSE
    Verbose level: from 0 (almost silent) to 5 (debug)
-p, --print 
    Print the parsed circuit
-o OUTFILE, --outfile=OUTFILE
    Data output file. Defaults to stdout.
--dc-guess=DC_GUESS 
    Guess to be used to start a OP or DC analysis: none or
    guess. Defaults to guess.
-t METHOD, --tran-method=METHOD
    Method to be used in transient analysis:
    implicit_euler, trap, gear2, gear3, gear4, gear5 or
    gear6. Defaults to TRAP.
--t-fixed-step 
    Disables the step control in transient analysis.
--v-absolute-tolerance=VEA
    Voltage absolute tolerance. Default: 1e-06 V
--v-relative-tolerance=VER
    Voltage relative tolerance. Default: 0.001
--i-absolute-tolerance=IEA
    Current absolute tolerance. Default: 1e-09 A
--i-relative-tolerance=IER
    Current relative tolerance. Default: 0.001
--h-min=HMIN 
    Minimum time step. Default: 1e-20
--dc-max-nr=DC_MAX_NR_ITER
    Maximum number of NR iterations for DC and OP
    analyses. Default: 10000
--t-max-nr=TRANSIENT_MAX_NR_ITER
    Maximum number of NR iterations for each time step
    during transient analysis. Default: 20
--t-max-time=TRANSIENT_MAX_TIME_ITER
    Maximum number of time iterations during transient
    analysis. Setting it to 0 (zero) disables the limit.
    Default: 0
--s-max-nr=SHOOTING_MAX_NR_ITER
    Maximum number of NR iterations during shooting
    analysis. Setting it to 0 (zero) disables the limit.
    Default: 10000
--gmin=GMIN 
    The minimum conductance to ground. Inserted when
    requested. Default: 1e-12
--cmin=CMIN 
    The minimum capacitance to ground. Default: 1e-18
