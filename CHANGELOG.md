<!-- Name: Changelog -->
# Changelog

## Sckzor's Changes

I have moderized the Python API calls to allow Akhab to continue to function
with modern python interpreters.

## Changelog for v 0.18

Version 0.18 represents the culmination of almost two months of efforts, yet
another step in the current time-based release scheme.

This release features support for the Fourier (.FOUR) and Fast Fourier Transform
directives (.FFT). Piece-wise linear functions are also now supported in Python,
with extra features like a repeat directory to enable periodicity.

Moreover, as in the previous releases, several commits were devoted to improving
the documentation: you can find the new, improved documentation online at
http://ahkab.readthedocs.org/en/latest/

### Features added:

* FEATURES

* Add .FOUR and .FFT support
* Add PWL waveforms (available from python only)

### Backwards incompatible changes:

* AC analyses are from now defined in Hz (compliant with SPICE). This also means
  that there is no omega (accessed through 'w') in the results set, it is now
  replaced by a frequency vector, expressed in Hz and accessed as 'f'.

User code can be usually easily fixed substituting ``2*pi*r[w]`` with ``r[f]``
or the like. Sorry for any trouble this creates.

### Changes from contributors and pull requests merged:

*None.* 

### Bugs fixed, short list:

BUGFIX: Fix opening of files in system that do not default to utf-8
BUGFIX: Suppress error message on Windows
BUGFIX: Use all assumptions, fix #32
BUGFIX: fix reading the netlist from stdin
BUGFIX: Use different options for SH and BFPSS
BUGFIX: no DC value -> return t=0
BUGFIX: convert to uppercase before comparison

## Changelog for v 0.17

Version 0.17 represents the culmination of 15 days of efforts, yet another step
in the current time-based release scheme.

This release features a new memoization system, AM and FM time functions,
improved docs and error message. We also got rid of the last instances of
`np.mat` in the code.

Moreover, as in the previous releases, several commits were devoted to improving
the documentation: you can find the new, improved documentation online at
http://ahkab.readthedocs.org/en/latest/

### Features added:

* FEATURES

* Implement a memoization system.
* Replace `mytime_function.value()` with `mytime_function(time)`.
* Add AM to the time functions.
* Add the SFFM time function.
* `find_vde_index()` now allows for `elem` and `part_id` args.
* `remove_elem()` now allows `elem` and `part_id` args.

### Backwards incompatible changes:

* the time functions have been moved to their own module, named
  `time_functions.py`. This means code referring to the `sin`, `pulse` and `exp`
  functions needs to change from `devices.sin` to `time_functions.sin`.
* You now need to call `set_results()` to set `pss_solution` data.

### Changes from contributors and pull requests merged:

*None.* 

### Bugs fixed, short list:

BUGFIX: plotting: remember to set the active figure.
BUGFIX: fix a bug in `circuit.remove_elem()`, add tests.
BUGFIX: fix a bug in `circuit.find_vde()`, add tests.
BUGFIX: fix error message in circuit.py.
BUGFIX: Do not ask for prediction if we don't have enough points.

## Changelog for v 0.16

Version 0.16 represents the culmination of 1 month of efforts, yet another step
in the current time-based release scheme.

### Features added:

This release brings the test coverage of the codebase over the 80% milestone for
the first time. It also features a faster substitution algorithm for symbolic
simulations, we empoly the ``tabulate`` module more extensively, producing a
prettier output (for example ``print_short()`` now prints a pretty table.)

Internally, We removed occurrences of ``np.matrix``, in favor of ``np.array``,
according to the overall planned switch in the library. This should be
completely transparent to the user.

A few changes that have been introduced are *incompatible* with the previous
releases:

* The axis iterators ``utilities.log_axis_iterator`` and
  ``utilities.lin_axis_iterator`` now follow the syntax ``(min, max, points)``.
* The solution method ``solution.asmatrix()`` has been renamed to
  ``solution.asarray()``.
* We droppped the (undocumented) support for accessing singularities in
  ``pz_solution`` as ``'Re(p0)'`` and ``'Im(p0)'``. Hopefully, being
  undocumented it had little use. Please use ``numpy.real(r['p0'])`` or
  ``numpy.imag(r['p0'])`` to achieve the same result.
* Remove ``printing.table_print()`` for ``print(printing.table())``.

We apologize about the above, we believe the technical debt we paid with the
changes above makes up for the discomfort to our userbase.

Moreover, as in the previous releases, several commits were devoted to improving
the documentation: you can find the new, improved documentation online at
http://ahkab.readthedocs.org/en/latest/

### Issues fixed

* #29 - ``ahkab`` should now work well in IPython running under Python2 again.

### Changes from contributors and pull requests merged:

*None.*

### Bugs fixed, short list:

* BUGFIX: Always plot in tests even if the test fails (especially then!)
* BUGFIX: catch ``ValueError`` in ``results.cid``
* BUGFIX: fix ``pss_solution.asmatrix()``
* BUGFIX: import ``codecs`` even if on Ipython
* BUGFIX: fix iterator off-by-one in ``symbolic_solution``
* BUGFIX: key misses in ``symbolic_solution`` raise ``KeyErrors``
* BUGFIX: fix slicing in the solution (use it too)
* BUGFIX: ``dc_solution.values()`` now slices along the correct axis
* BUGFIX: ``load_csv()`` raises ``ValueError``
* BUGFIX: key misses raise ``KeyErrors``
* BUGFIX: prevent breaking readline on ipython, fixes #29
* BUGFIX: ensure w is always returned as real data, not cplx
* BUGFIX: Improve the iterators syntax, always include endpoints.
* BUGFIX: raise ``KeyError`` for key misses in ac_solution
* BUGFIX: ``items()`` returns no arrays.
* BUGFIX: ``values()`` returns a list, not an array.
* BUGFIX: fix iterator increment-by-one bug

## Changelog for v 0.15

Version 0.15 is just a bugfix release addressing the fact that the minimum version of the ``sympy`` release that is needed is 0.7.6.

## Changelog for v 0.14

Version 0.14 represents the culmination of three months of efforts, yet another step in the current time-based release scheme.

This release features much improved tabular print-outs, thanks to ``tabulate``, faster circuit checking and reverse node look-ups.

Moreover, as in the previous releases, many commits were devoted to improving the documentation: you can find the new, improved documentation online at http://ahkab.readthedocs.org/en/latest/

### Features added:

* Use ``tabulate`` to print out pretty tables. *They are pretty, oh so pretty tables.*
* Add and use ``circuit.get_nodes_number()``.
* Change ``nodes_dict`` to speed up reverse lookups.
* Faster duplicate IDs check.

### Changes from contributors and pull requests merged:

*None.*

### Bugs fixed, short list:

* BUGFIX: row/col reference,
* BUGFIX: Suppress printing of the netlist file if there is None.
* BUGFIX: Add a wrapper around stdout to fix encoding errors on *UNIX, when no locale is set (PY2 only).
* BUGFIX: Print warning if the locale is not set.


## Changelog for v 0.13

Version 0.13 represents the culmination of efforts dating back to January and contains 135 commits.

This release features current-controlled current sources and voltage sources and new options to select the format of numbers printed for user display and how the MNA matrix is built in Symbolic Analysis.

Moreover, as in the previous releases, many commits were devoted to improving the documentation: you can find the new, improved documentation at http://ahkab.readthedocs.org/en/latest/

### Features added:

* Implement CCCS, fixes #22
* Implement CCVS, fixes #21
* Add option to formulate the MNA and N matrices with resistors for symbolic analysis.
* Add options to suppress zeros and select printing precision.

### Changes from contributors and pull requests merged:

*None.*

### Bugs fixed, short list:

* BUGFIX: fix printing of missing elements.
* BUGFIX: fix missing import
* BUGFIX: reset look-up table when the temperature changes
* BUGFIX: catch missing matplotlib
* BUGFIX: remove early test It masks more detailed tests below.
* BUGFIX: fix get_netlist_elem_line()Pass the element nodes through nodes_dict
* BUGFIX: correct the netlist label for tran start time (tstart)
* BUGFIX: fix parsing of capacitors wrt comments.
* BUGFIX: fix all Symbols to be uppercase
* BUGFIX power by VCCS
* BUGFIX: Fix total power computation for GIsources.

## Changelog for v 0.12

Version 0.12 represents the culmination of efforts dating back to May and contains 142 commits.

This release has no new features, but it has several bug fixes and it introduces Python3 support (!).

Moreover, many commits were devoted to improving the documentation: you can find the new, improved documentation at http://ahkab.readthedocs.org/en/latest/

### Features added:

*None.*

### Changes from contributors and pull requests merged:

*None.*

### Bugs fixed, short list:

* BUGFIX: Respect the user config when saving to file.
* BUGFIX: Respect CLI-specified transient DFs.
* BUGFIX: use the warning interface.
* BUGFIX: always set the simulation options before running tests.
* BUGFIX: fix detection of wd
* BUGFIX: .include now recovers the path relative to the netlist
* BUGFIX: fix printing of .DC statements.
* BUGFIX: fix comparison of int and NoneType (new Python behaviour)


## Changelog for v 0.11

Version 0.11 represent the culmination of four months of efforts.

### Features added:

* Add Pole-Zero (.pz) analysis.
* Add sparse matrix support and LU factorization with SuperLU.
* AC solutions data is returned as complex floats.
* The brute-force PSS algo is now known as 'brute-force'.

Ahkab also got considerably faster in DC/OP/TRAN/PSS analyses.

A PY3 version passing the whole test suite is available in the
repository.

### Changes from contributors and pull requests merged:

*None.*

### Bugs fixed, short list:

* BUGFIX: AC sweeps consistently default to LOG.
* BUGFIX: Convert print_netlist_elem_line to get_...
* BUGFIX: parse_analysis() is not a generator anymore`.
* BUGFIX: Remove print_circuit(), use print(circuit).
* BUGFIX: EKV devices now match the usual constructor.
* BUGFIX: do not attempt to print conv. details if sing. MNA
* BUGFIX: Use ndarrays (almost) everywhere, drop mats.
* BUGFIX: Don't crash if a plot label is not found.
* BUGFIX: CMIN-related crash.
* BUGFIX: Remove hack for now closed bug in sympy.
* BUGFIX: Reset sim options before skipping a test.

## Changelog for v 0.10

### Bugs fixed, short list:

*    BUGFIX: allow printing of PSS results.
*    BUGFIX: inductor coupling and time function parsing.
*    BUGFIX: fix checking of symbolic results.
*    BUGFIX: InductorCoupling: init syntax and `add print_netlist_elem_line()`.
*    BUGFIX: fix regression in TF calculation (wrong sympy assumptions).
*    BUGFIX: function parsing
*    BUGFIX: fix 'autonomous' dest value.
*    BUGFIX: `modify_x0_for_ic()` should NOT modify op results.
*    BUGFIX: op table size
*    BUGFIX: MOSQ init param order
*    BUGFIX: Symbolic analysis bugfix.
*    BUGFIX: Fix DC sweeps.
*    BUGFIX: Fix dict access of OP results. OP info now goes to info file.
*    BUGFIX: Deep copy `an_list` so that we can pop out items as we please transparently to the user.
*    BUGFIX: take the file.op into account. Also abs paths.
*    BUGFIX: load_csv crash, add late check for headers consistency.
*    BUGFIX: Fix the sine function. Notice the meaning of theta changed with this commit.
*    BUGFIX: Correct assumptions in subs.
*    BUGFIX: fix `element.part_id` in subckts
*    BUGFIX: testing, do not compare time exec results from different machines.

### Changes from contributors and pull requests merged:

*    Incorporate most changes from @weilawei.
*    Merge pull request #18 from endolith/patch-1
*    Merge branch 'patch-1' of https://github.com/mightyiam/ahkab

### Features added:

*    Use `__version__.py` for versioning
*    Diode: Enable caching. Simplify interfaces.
*    Allow running lists of analyses, or 1 analysis alone too.
*    Add gnd as `mycirc.gnd`.
*    Add command line entry point to `setup.py`.
*    Testing: Add APITest to automate testing of the `ahkab` Python interface.
*    Testing: Test suite running on Travis, with coverage.
*    Testing: Selective test skipping on Travis-CI

### Changes:

*    Massive refactoring of the code base: enforce PEP8 everywhere.
*    Symbolic: transfer functions are results instances too.
*    Plotting: remove useless functions, keep some as internal-only.
*    `get_headers_index(..., load_headers=[])` returns all header indices.
*    Testing: No reference runs with nose.
*    Ticker: simplify the ticker API.
*    Speed up the matrix assembly precomputing conductances
*    Allow printing out to stdout also `if filename == '-' or sys.stdout`
*    Matplotlib auto backend selection for machines without `$DISPLAY`
*    Rewrite csvlib to use numpy's R&W methods
*    Elements now print out their own netlist entries.
*    Move `convergence_check()` to utilities.
*    Use relative imports and move to single versioning of files.
*    Remove def. values from time functions

Thanks to @weilawei**, @endolith, @mightyiam for their contributions.

** Big thanks man, excellent job.

## Changes in previous versions, in brief:

## 2/0913
*  Added support for voltage-controlled switches

## 20/08/13
*  Improved diode model with temperature support (but no capacitances).

## 8/8/13
 * Dropped psycho support.
 * Re-added the quadratic law MOS transistor model.
 * Updated wiki pages and examples.
 * Many small fixes.

## 20/12/11
* Ahkab can now be run in a Python program.

### 04/01/11
* Dropped gnuplot support, ahkab now requires matplotlib.
