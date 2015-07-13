Ahkab |release| Docs
====================

:Release: |release|
:Date: |today|

Ahkab (pronounced "uh, cab") is an open-source SPICE-like interactive circuit simulator.

|Build Status| |Coverage Status| |PyPi version| |GPLv2 license| |DOI|

.. these entries are hidde, they go to the sidebar but they do not clutter the
   main page.

.. toctree::
   :maxdepth: 1
   :hidden:

   Homepage <https://ahkab.github.io/ahkab/>
   Repository <https://github.com/ahkab/ahkab>
   Issue Tracker <https://github.com/ahkab/ahkab/issues>


General help pages
------------------

.. toctree::
   :maxdepth: 1

   misc/A-Good-Idea
   help/Install-Notes
   help/Command-Line-Help
   help/Netlist-Syntax

Getting started: examples and tutorials
---------------------------------------

``ahkab`` can be used from Python as a module and from the shell through its
Command Line Interface (CLI).


Simulating from Python
**********************

.. toctree::
   :maxdepth: 1

   examples/Simple_OP
   examples/Python_API
   examples/PZ_Example

Simulating from the command line
********************************

.. toctree::
   :maxdepth: 1

   examples/OP_simulation
   examples/Transient-Example
   examples/Symbolic-simulation

Help pages on particular elements
---------------------------------

.. toctree::
   :maxdepth: 1

   help/Mutual-Inductors

Module reference
----------------

.. toctree::
   :maxdepth: 1

   ahkab
   ac
   bfpss
   circuit
   constants
   csvlib
   dc_analysis
   dc_guess
   devices
   diode
   ekv
   fourier
   gear
   implicit_euler
   mosq
   netlist_parser
   options
   plotting
   printing
   pss
   pz
   results
   shooting
   switch
   symbolic
   testing
   ticker
   time_functions
   transient
   trap
   utilities

License
-------

.. toctree::
   :maxdepth: 1

   misc/COPYING

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

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
