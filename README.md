#<img src="http://raw.github.com/wiki/ahkab/ahkab/images/logo_small.png" alt="Monkeying around" style="width: 80px;"/> ahkab

***a SPICE-like electronic circuit simulator written in Python***

The code should be easy to read and modify, the main language is Python 2.x and it is platform-independent.[![githalytics.com alpha](https://cruel-carlota.pagodabox.com/3f4b146d6a15f66802f1906e5cf4f68c "githalytics.com")](http://githalytics.com/ahkab/ahkab)

## News! ##

 * Ahkab v0.11 was released on Oct 22 2014, including several bugfixes and improvements. It is recommended to upgrade. Check out [the release notes](https://github.com/ahkab/ahkab/releases/tag/v0.11) for more!
 * The whole codebase has been going through a (yet incomplete) refactoring and documenting effort. The [new documentation is available on RTD](http://ahkab.readthedocs.org/en/latest/).
 * Tests have been added to the code base, to weed out bugs and have a consistent way to check for regressions. If you have tests you would like to suggest, please send a PR my way.
 *  My resources are limited these days, so the much-needed work is proceeding slowly, albeit hopefully steadily. If you are interested and you would like to contribute to refactoring or documenting a particular feature, it would be very welcome.

[![Build Status](https://travis-ci.org/ahkab/ahkab.png?branch=master)](https://travis-ci.org/ahkab/ahkab) [![Coverage Status](https://coveralls.io/repos/ahkab/ahkab/badge.png?branch=master)](https://coveralls.io/r/ahkab/ahkab?branch=master)
[![PyPi version](http://img.shields.io/badge/version-0.11-brightgreen.png)](https://pypi.python.org/pypi/ahkab/) [![GPLv2 license](http://img.shields.io/badge/license-GPL%20v2-brightgreen.png)](https://raw.githubusercontent.com/ahkab/ahkab/master/LICENSE)
<!--- [![PyPi downloads](https://pypip.in/download/ahkab/badge.png)](https://pypi.python.org/pypi/ahkab/) --->


## Supported simulations: ##
  * Numeric:
    * **Operating point**, with guess computation to speed up the solution. See example: [Downscaling current mirror](https://github.com/ahkab/ahkab/wiki/Example:-OP-simulation)
    * **DC sweep**
    * **Transient analysis**, available differentiation formulas: implicit Euler, trapezoidal, gear orders from 2 to 5. See for example the [simulation of a Colpitts Oscillator](https://github.com/ahkab/ahkab/wiki/Example:-Transient-simulation-1).
    * **AC analysis**
    * **PZ** analysis
    * **Periodic steady state analysis** of non-autonomous circuits, _time_ _domain_ shooting and brute-force algorithms.
  * Symbolic: 
    * **Small signal analysis**, AC or DC, with extraction of transfer functions, DC gain, poles and zeros. Various [symbolic analysis examples on this page](https://github.com/ahkab/ahkab/wiki/Example:-Symbolic-simulation).

The results are saved to disk, plotted or printed to stdout and can be read/processed by the most common tools (eg. [Octave](http://www.gnu.org/software/octave/), [gnuplot](http://www.gnuplot.info/), [Matlab](http://www.mathworks.com/products/matlab/), [gwave](http://www.telltronics.org/software/gwave/) and others)

##Install##

The program requires:
* the **Python 2 interpreter** (at least v.2.6, and Python3 is in the works),
* numpy>=1.7.0, matplotlib and sympy.

If you need more information about the dependencies, check the [Install notes](https://github.com/ahkab/ahkab/wiki/Install:-Notes).

##Usage##

 * `ahkab` can be run _within Python scripts as a library_. 

###Example###
 
<img src="https://rawgithub.com/ahkab/ahkab/master/doc/images/readme_example/pbf.svg" alt="Example schematic: a 5th order 1kHz band-pass Butterworth filter"/>

```python
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
```

<img src="https://rawgithub.com/ahkab/ahkab/master/doc/images/readme_example/bpf_results.svg" alt="Example: AC simulation resultsr"/>

 * or stand-alone with a netlist file, the syntax being:

    `$ python ahkab -o graph.dat <netlist file>`

See `ahkab --help` for command line switches.

###Documentation###

Refer to the [netlist syntax page](https://github.com/ahkab/ahkab/wiki/Help:-Netlist-Syntax) for how to write the netlist files that describe the circuit. Experience with running SPICE or other commercial simulators can be useful.

The latter option is shown briefly in the **[Python API](https://github.com/ahkab/ahkab/wiki/Example:-Python-API)** wiki page. The code comes with docstrings associated with _most_ functions, type `help(ahkab.function_name)`.

### How this project was born ###

This project was born when I was an enthusistic undergrad, apparently with plenty of free time, attending "Simulazione Circuitale" (*Circuit Simulation*) taught by [Prof. A. Brambilla](http://brambilla.dei.polimi.it/) back in Italy at the Polytechnic University of Milan.

I am grateful to prof. Brambilla for teaching one of the most interesting courses of my university years. -GV

### Contributors ###
[Giuseppe Venturini](https://github.com/ggventurini), [Ian Daniher](https://github.com/itdaniher), [Rob Crowther](https://github.com/weilawei).

### Bugs and patches ###

Note that _I often add new functionality at the expense of breaking stuff_. Most likely I will introduce a new feature even if that means breaking a couple of others. It should get fixed soon, but if you have a bit of time to spare, you can send me a pull request or a patch. :)

Does it work? Bugs? Do you have patches? Did you run some noteworthy simulation? Let me know! Feedback is very welcome, my [email address](http://tinymailto.com/5310) is available after a captcha.
