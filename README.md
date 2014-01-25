#<img src="http://raw.github.com/wiki/ahkab/ahkab/images/logo_small.png" alt="Monkeying around" style="width: 80px;"/> ahkab

***a SPICE-like electronic circuit simulator written in Python***

The code should be easy to read and modify, the main language is Python 2.x and it is platform-independent.[![githalytics.com alpha](https://cruel-carlota.pagodabox.com/3f4b146d6a15f66802f1906e5cf4f68c "githalytics.com")](http://githalytics.com/ahkab/ahkab)

### News! ###

 * **Update your scientific Python installation!** We benefit from the latest updates to numpy and get **a 40% speed-up** for all large networks simulations.

![r2r](https://f.cloud.github.com/assets/5038142/1903519/85895ee4-7c83-11e3-9712-c63c7152ce23.png)


### Supported simulations: ###
  * Numeric:
    * **Operating point**, with guess computation to speed up the solution. See example: [Downscaling current mirror](https://github.com/ahkab/ahkab/wiki/Example:-OP-simulation)
    * **DC sweep**
    * **Transient analysis**, available differentiation formulas: implicit Euler, trapezoidal, gear orders from 2 to 5. See for example the [simulation of a Colpitts Oscillator](https://github.com/ahkab/ahkab/wiki/Example:-Transient-simulation-1).
    * **AC analysis**
    * **Periodic steady state analysis** of non-autonomous circuits, _time_ _domain_ shooting and brute-force algorithms.
  * Symbolic: 
    * **Small signal analysis**, AC or DC, with extraction of transfer functions, DC gain, poles and zeros. Various [symbolic analysis examples on this page](https://github.com/ahkab/ahkab/wiki/Example:-Symbolic-simulation).

The results are saved to disk, plotted or printed to stdout and can be read/processed by the most common tools (eg. [Octave](http://www.gnu.org/software/octave/), [gnuplot](http://www.gnuplot.info/), [Matlab](http://www.mathworks.com/products/matlab/), [gwave](http://www.telltronics.org/software/gwave/) and others)

###Download and install###

There are no packages for the time being (this program is at an early development stage). Go to [ahkab on github](https://github.com/ahkab/ahkab) and follow the instructions to check out the code. You can find the list of the dependencies in the [Install notes](https://github.com/ahkab/ahkab/wiki/Install:-Notes).

###Usage###

 * `ahkab` can be run _within Python scripts as a library_. This will likely become the preferred way in the future. See **[this butterworth filter simulation](https://github.com/ahkab/ahkab/wiki/Example:-Python-API)** for an example/tentative tutorial.

 * or stand-alone with a netlist file, the syntax being:

    $ python ahkab -o graph.dat <netlist file>`

See `ahkab --help` for command line switches.

###Documentation###

The simulator can either be run from the command line with a netlist file or included in a python script. Both possibilities will be maintained for the foreseeable future. 

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
