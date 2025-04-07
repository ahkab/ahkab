ahkab
=====
Updated to be compatible with python3.12 and beyond, and scipy 1.x
=====

**a SPICE-like electronic circuit simulator written in Python**

The code should be easy to read and modify, the main language is Python
-- 2 or 3 -- and it is platform-independent.

News!
-----

-  Ahkab v0.18 was released on July 18 2015, including new features,
   bugfixes and improved documentation. It is recommended to upgrade.
   Check out `the release
   notes <https://github.com/ahkab/ahkab/releases/tag/v0.18>`__ for
   more!
-  The whole codebase has been going through a (yet incomplete)
   refactoring and documenting effort. The `new documentation is
   available on RTD <http://ahkab.readthedocs.org/en/latest/>`__.

My resources are limited these days, so the much-needed work is
proceeding slowly, albeit hopefully steadily. If you are interested and
you would like to contribute to refactoring or documenting a particular
feature, it would be very welcome.

|Build Status| |Coverage Status| |PyPi version| |GPLv2 license| |DOI|

Supported simulations:
----------------------

-  Numeric:

   -  **Operating point**, with guess computation to speed up the
      solution. See example: `Downscaling current
      mirror <https://ahkab.readthedocs.org/en/latest/examples/OP_simulation.html>`__
   -  **DC sweep**
   -  **Transient analysis**, available differentiation formulas:
      implicit Euler, trapezoidal, gear orders from 2 to 5. See for
      example the `simulation of a Colpitts
      Oscillator <https://ahkab.readthedocs.org/en/latest/examples/Transient-Example.html>`__.
   -  **AC analysis**
   -  **PZ** analysis
   -  **Periodic steady state analysis** of non-autonomous circuits,
      *time* *domain* shooting and brute-force algorithms.

-  Symbolic:

   -  **Small signal analysis**, AC or DC, with extraction of transfer
      functions, DC gain, poles and zeros. Various `symbolic analysis
      examples on this
      page <https://ahkab.readthedocs.org/en/latest/examples/Symbolic-simulation.html>`__.

The results are saved to disk, plotted or printed to stdout and can be
read/processed by the most common tools (eg.
`Octave <http://www.gnu.org/software/octave/>`__,
`gnuplot <http://www.gnuplot.info/>`__,
`Matlab <http://www.mathworks.com/products/matlab/>`__,
`gwave <http://www.telltronics.org/software/gwave/>`__ and others)

Install
-------

The program requires:

-  the Python interpreter version 2 or 3 (at least v.2.6 for Python2,
   v.3.3 for Python3),
-  numpy>=1.7.0, scipy>=0.14.0, sympy>=0.7.6 and tabulate>=0.7.3.

Matplotlib is strongly recommended and no plotting will work without.

If you need more information about the dependencies, check the `Install
notes <https://ahkab.readthedocs.org/en/latest/help/Install-Notes.html>`__.

Usage
-----

1. ``ahkab`` can be run as a Python library
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

2. ``ahkab`` can be run from the command line with a netlist file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The syntax is:

::

    `$ python ahkab -o graph.dat <netlist file>`

See ``ahkab --help`` for command line switches, `also online on the
documentation
pages. <http://ahkab.readthedocs.org/en/latest/help/Command-Line-Help.html>`__

Documentation
~~~~~~~~~~~~~

The `documentation is available on
RTD <http://ahkab.readthedocs.org/en/latest/>`__.

There, you can find a
`documentation <http://ahkab.readthedocs.org/en/latest/ahkab.html>`__
and
`examples <http://ahkab.readthedocs.org/en/latest/examples/Python_API.html>`__
regarding how to simulate from a Python script.

Refer to the `netlist syntax
page <http://ahkab.readthedocs.org/en/latest/help/Netlist-Syntax.html>`__
if you prefer to write netlist files that describe the circuit.

Experience with running SPICE or related commercial simulators can be
very useful: this is not for the faint of heart.

Development model
~~~~~~~~~~~~~~~~~

-  The development happens on the `github
   repository <https://github.com/ahkab/ahkab>`__,
-  Mostly on the master branch, with feature branch being created only
   for special purposes or non-trivial features.
-  Snapshots are released on a (hopefully) regular basis and are
   available on the `Releases pages, complete with
   changelog <https://github.com/ahkab/ahkab/releases>`__ and on
   `PYPI <https://pypi.python.org/pypi/ahkab/>`__

Patches and pull requests are welcome!

How this project was born
~~~~~~~~~~~~~~~~~~~~~~~~~

This project was born when I was an enthusistic undergrad, apparently
with plenty of free time, attending "Simulazione Circuitale" (*Circuit
Simulation*) taught by `Prof. A.
Brambilla <http://brambilla.dei.polimi.it/>`__ back in Italy at the
Polytechnic University of Milan.

I am grateful to prof. Brambilla for teaching one of the most
interesting courses of my university years. -GV

Bugs and patches
~~~~~~~~~~~~~~~~

Does it work? Bugs? Do you have patches? Did you run some noteworthy
simulation? Let me know! Feedback is very welcome, my `email
address <http://tinymailto.com/5310>`__ is available after a captcha.

Support the development with a donation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you wish to support the development of ``ahkab``, ***please donate to
cancer research:***

-  `Association for International Cancer Research
   (eng) <http://www.aicr.org.uk/donate.aspx>`__,

or

-  `Fond. IRCCS Istituto Nazionale dei Tumori
   (it) <http://www.istitutotumori.mi.it/modules.php?name=Content&pa=showpage&pid=24>`__.

Credits
~~~~~~~

**Authors:** `Giuseppe Venturini <https://github.com/ggventurini>`__,
with contributions from `Ian Daniher <https://github.com/itdaniher>`__
and `Rob Crowther <https://github.com/weilawei>`__.

**Code:** the module ``py3compat.py`` is (c) 2013 - the Jinja team.

**Dependencies:** many thanks to the authors of ``numpy``, ``scipy``,
``sympy``, ``matplotlib`` and ``tabulate``!

.. |Build Status| image:: https://travis-ci.org/ahkab/ahkab.png?branch=master
   :target: https://travis-ci.org/ahkab/ahkab
.. |Coverage Status| image:: https://coveralls.io/repos/ahkab/ahkab/badge.png?branch=master
   :target: https://coveralls.io/r/ahkab/ahkab?branch=master
.. |PyPi version| image:: http://img.shields.io/badge/version-0.18-brightgreen.png
   :target: https://pypi.python.org/pypi/ahkab/
.. |GPLv2 license| image:: http://img.shields.io/badge/license-GPL%20v2-brightgreen.png
   :target: https://raw.githubusercontent.com/ahkab/ahkab/master/LICENSE
.. |DOI| image:: https://zenodo.org/badge/doi/10.5281/zenodo.19967.svg
   :target: http://dx.doi.org/10.5281/zenodo.19967
