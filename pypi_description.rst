 ahkab
======

***a SPICE-like electronic circuit simulator written in Python***

The code should be easy to read and modify, the main language is Python
-- 2 or 3 -- and it is platform-independent.\ |githalytics.com alpha|

News!
-----

-  Happy new year! 2015 was kicked off with a brand new Ahkab release.
   Ahkab v0.12 was released on Jan 05 2015, including several bugfixes,
   improvements, a revised documentation and Python3 support. It is
   recommended to upgrade. Check out `the release
   notes <https://github.com/ahkab/ahkab/releases/tag/v0.12>`__ for
   more!
-  The whole codebase has been going through a (yet incomplete)
   refactoring and documenting effort. The `new documentation is
   available on RTD <http://ahkab.readthedocs.org/en/latest/>`__.

My resources are limited these days, so the much-needed work is
proceeding slowly, albeit hopefully steadily. If you are interested and
you would like to contribute to refactoring or documenting a particular
feature, it would be very welcome.

|Build Status| |Coverage Status| |PyPi version| |GPLv2 license|

Supported simulations:
----------------------

-  Numeric:

   -  **Operating point**, with guess computation to speed up the
      solution. See example: `Downscaling current
      mirror <https://github.com/ahkab/ahkab/wiki/Example:-OP-simulation>`__
   -  **DC sweep**
   -  **Transient analysis**, available differentiation formulas:
      implicit Euler, trapezoidal, gear orders from 2 to 5. See for
      example the `simulation of a Colpitts
      Oscillator <https://github.com/ahkab/ahkab/wiki/Example:-Transient-simulation-1>`__.
   -  **AC analysis**
   -  **PZ** analysis
   -  **Periodic steady state analysis** of non-autonomous circuits,
      *time* *domain* shooting and brute-force algorithms.

-  Symbolic:

   -  **Small signal analysis**, AC or DC, with extraction of transfer
      functions, DC gain, poles and zeros. Various `symbolic analysis
      examples on this
      page <https://github.com/ahkab/ahkab/wiki/Example:-Symbolic-simulation>`__.

The results are saved to disk, plotted or printed to stdout and can be
read/processed by the most common tools (eg.
`Octave <http://www.gnu.org/software/octave/>`__,
`gnuplot <http://www.gnuplot.info/>`__,
`Matlab <http://www.mathworks.com/products/matlab/>`__,
`gwave <http://www.telltronics.org/software/gwave/>`__ and others)

Install
-------

The program requires: \* the Python interpreter version 2 or 3 (at least
v.2.6 for Python2, v.3.3 for Python3), \* numpy>=1.7.0, scipy>0.11.1,
matplotlib and sympy.

If you need more information about the dependencies, check the `Install
notes <https://github.com/ahkab/ahkab/wiki/Install:-Notes>`__.

Usage
-----

-  ``ahkab`` can be run *within Python scripts as a library*.

Example
~~~~~~~

.. code:: python

    from ahkab import new_ac, run
    from ahkab.circuit import Circuit
    from ahkab.plotting import plot_results # calls matplotlib for you
    import numpy as np

    # Define the circuit
    cir = Circuit('Butterworth 1kHz band-pass filter')
    cir.add_vsource('V1', 'n1', cir.gnd, dc_value=0., ac_value=1.)
    cir.add_resistor('R1', 'n1', 'n2', 50.)
    cir.add_inductor('L1', 'n2', 'n3', 0.245894)
    cir.add_capacitor('C1', 'n3', 'n4', 1.03013e-07)
    cir.add_inductor('L2', 'n4', cir.gnd, 9.83652e-05)
    cir.add_capacitor('C2', 'n4', cir.gnd, 0.000257513)
    cir.add_inductor('L3', 'n4', 'n5', 0.795775)
    cir.add_capacitor('C3', 'n5', 'n6', 3.1831e-08)
    cir.add_inductor('L4', 'n6', cir.gnd, 9.83652e-05)
    cir.add_capacitor('C4', 'n6', cir.gnd, 0.000257513)
    cir.add_capacitor('C5', 'n7', 'n8', 1.03013e-07)
    cir.add_inductor('L5', 'n6', 'n7', 0.245894)
    cir.add_resistor('R2', 'n8', cir.gnd, 50.)

    # Define the analysis
    ac1 = new_ac(2.*np.pi*.97e3, 2.*np.pi*1.03e3, 1e2, x0=None)

    # run it
    res = run(cir, ac1)

    # plot the results
    plot_results('5th order 1kHz Butterworth filter', [('|Vn8|',"")], res['ac'],
                 outfilename='bpf_transfer_fn.png')

-  or stand-alone with a netlist file, the syntax being:

   ``$ python ahkab -o graph.dat <netlist file>``

See ``ahkab --help`` for command line switches.

Documentation
~~~~~~~~~~~~~

The `new documentation is available on
RTD <http://ahkab.readthedocs.org/en/latest/>`__.

There, you can find a
`documentation <http://ahkab.readthedocs.org/en/latest/ahkab.html>`__
and
`examples <http://ahkab.readthedocs.org/en/latest/examples/Python_API.html>`__
regarding how to simulate from a Python script.

Refer to the `netlist syntax
page <http://ahkab.readthedocs.org/en/latest/help/Netlist-Syntax.html>`__
if you prefer to write netlist files that describe the circuit.

Experience with running SPICE or other commercial simulators can be very
useful: this is not for the faint of heart.

How this project was born
~~~~~~~~~~~~~~~~~~~~~~~~~

This project was born when I was an enthusistic undergrad, apparently
with plenty of free time, attending "Simulazione Circuitale" (*Circuit
Simulation*) taught by `Prof. A.
Brambilla <http://brambilla.dei.polimi.it/>`__ back in Italy at the
Polytechnic University of Milan.

I am grateful to prof. Brambilla for teaching one of the most
interesting courses of my university years. -GV

Credits
~~~~~~~

**Contributors:** `Giuseppe
Venturini <https://github.com/ggventurini>`__, `Ian
Daniher <https://github.com/itdaniher>`__, `Rob
Crowther <https://github.com/weilawei>`__.

**Code:** the module ``py3compat.py`` is (c) 2013 - the Jinja team.

**Deps:** many thanks to the authors of ``numpy``, ``scipy`` and
``sympy``!

Bugs and patches
~~~~~~~~~~~~~~~~

Note that *I often add new functionality at the expense of breaking
stuff*. Most likely I will introduce a new feature even if that means
breaking a couple of others. It should get fixed soon, but if you have a
bit of time to spare, you can send me a pull request or a patch. :)

Does it work? Bugs? Do you have patches? Did you run some noteworthy
simulation? Let me know! Feedback is very welcome, my `email
address <http://tinymailto.com/5310>`__ is available after a captcha.

.. |githalytics.com alpha| image:: https://cruel-carlota.pagodabox.com/3f4b146d6a15f66802f1906e5cf4f68c
   :target: http://githalytics.com/ahkab/ahkab
.. |Build Status| image:: https://travis-ci.org/ahkab/ahkab.png?branch=master
   :target: https://travis-ci.org/ahkab/ahkab
.. |Coverage Status| image:: https://coveralls.io/repos/ahkab/ahkab/badge.png?branch=master
   :target: https://coveralls.io/r/ahkab/ahkab?branch=master
.. |PyPi version| image:: http://img.shields.io/badge/version-0.12-brightgreen.png
   :target: https://pypi.python.org/pypi/ahkab/
.. |GPLv2 license| image:: http://img.shields.io/badge/license-GPL%20v2-brightgreen.png
   :target: https://raw.githubusercontent.com/ahkab/ahkab/master/LICENSE
