<!-- Name: Changelog -->
#Changelog

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
