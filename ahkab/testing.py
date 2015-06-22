# -*- coding: utf-8 -*-
# testing.py
# Testing framework
# Copyright 2014 Giuseppe Venturini

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

"""
A straight-forward framework to buid tests to ensure no regressions
occur during development.

Two classes for describing tests are defined in this module:

- :class:`NetlistTest`, used to run a netlist-based test,
- :class:`APITest`, used to run an API-based test.

Every test, no matter which class is referenced internally, is
univocally identified by a alphanumeric id, which will
be referred to as ``<test_id>`` in the following.


Directory structure
\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"

The tests are placed in ``tests/``, under a directory with the same
id as the test, ie:

::

    tests/<test_id>/


Running tests
\"\"\"\"\"\"\"\"\"\"\"\"\"

The test is performed with as working directory one among the following:

 - The ahkab repository root,

 - ``tests/``,

 - ``tests/<test_id>``.

this is necessary for the framework to find its way to the reference files.

More specifically a test can either be run manually through the Python
interpreter:

::

    python tests/<test_id>/test_<test_id>.py

or with the ``nose`` testing package:

::

    nosetests tests/<test_id>/test_<test_id>.py

To run the whole test suite, issue:

::

    nosetests tests/*/*.py

Please refer to the `nose documentation`_ for more info about the command
``nosetests``.

.. _nose documentation: https://nose.readthedocs.org/en/latest/

Running your tests for the first time
\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"

The first time you run a test you defined yourself, no reference data will be
available to check the test results and decide whether the test was passed or
if a test fail occurred.

In this case, if you call ``nose``, the test will (expectedly) fail.

Please run the test manually (see above) and the test framework will generate
the reference data for you.

Please *check the generated reference data carefully!*
Wrong reference defeats the whole concept of running tests!


Overview of a typical test based on :class:`NetlistTest`
\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"

Each test is composed by multiple files.

Required files
^^^^^^^^^^^^^^

The main directory must contain:

- ``<test_id>.ini``, an INI configuration file containing the details of the
  test,

- ``test_<test_id>.py``, the script executing the test,

- ``<test_id>.ckt``, the main netlist file to be run.

- the reference data files for checking the pass/fail status of the test.
  These can be automatically generated, as it will be shown below.

With the exception of the netlist file, which is free for the test writer
to define, and the data files, which clearly depend on the test at hand,
the other files have a predefined structure which will be examined
in more detail in the next sections.

Configuration file
''''''''''''''''''

Few rules are there regarding the entries in the configuration file.

They are as follows:

- The file name must be ``<test_id>.ini``,

- It must be located under ``tests/<test_id>/``,

- It must have a ``[test]`` section, containing the following entries:

  - ``name``, set to the ``<test_id>``, for error-checking,

  - ``netlist``, set to the netlist filename, ``<test_id>.ckt``, prepended
    with the the netlist path relative to ``tests/<test_id>/`` (most of
    the time that means just ``<test_id>.ckt``)

  - ``type``, a comma-separated list of analyses that will be executed during
    the test. Values may be ``op``, ``dc``, ``tran``, ``symbolic``... and so on.

  - One entry ``<analysis>_ref`` for each of the analyses listed in the
    ``type`` entry above.
    The value is recommended to be set to ``<test_id>-ref.<analysis>`` or
    ``<test_id>-ref.<analysis>.pickle``, if you prefer to save data in
    Python's pickle format. Notice only trusted pickle files should
    ever be loaded.

  - ``skip-on-travis``, set to either ``0`` or ``1``, to flag whether this
    test should be run on Travis-CI or not. Torture tests, tests needing
    lots of CPU or memory, and long-lasting tests in general should be
    disabled on Travis-CI to not exceed:

    - a total build time of 50 minutes,

    - A no stdout activity time of 10 minutes.

  - ``skip-on-pypy``, set to either ``0`` or ``1``, to flag whether the test
    should be skipped if useing a PYPY Python implemetntation or not. In
    general, as PYPY supports neither ``scipy`` nor ``matplotlib``, only
    symbolic-oriented tests make sense with PYPY (where it really excels!).


The contents of an example test configuration file ``rtest1.ini``
follow, as an example.

::

    [test]
    name = rtest1
    netlist = rtest1.ckt
    type = dc, op
    dc_ref = rtest1-ref.dc
    op_ref = rtest1-ref.op
    skip-on-travis = 0
    skip-on-pypy = 1


Script file
'''''''''''

The test script file is where most of the action takes place and where
the highest amount of flexibility is available.

That said, the ahkab testing framework was designed to make for extremely
simple and straight-forward test scripts.

It is probably easier to introduce writing the scripts with an example.

Below is a typical script file.

::

    from ahkab.testing import NetlistTest
    from ahkab import options
    # add this to prevent interactive plot directives
    # in the netlist from halting the test waiting for
    # user input
    options.plotting_show_plots = False

    def myoptions():
        # optionally, set non-standard options
        sim_opts = {}
        sim_opts.update({'gmin':1e-9})
        sim_opts.update({'iea':1e-3})
        sim_opts.update({'transient_max_nr_iter':200})
        return sim_opts

    def test():
        # this requires a netlist ``mytest.ckt``
        # and a configuration file ``mytest.ini``
        nt = NetlistTest('mytest', sim_opts=myoptions())
        nt.setUp()
        nt.test()
        nt.tearDown()

    # It is recommended to set the docstring to a meaningful value
    test.__doc__ = "My test description, printed out by nose"

    if __name__ == '__main__':
        nt = NetlistTest('mytest', sim_opts=myoptions())
        nt.setUp()
        nt.test()

Notice how a function ``test()`` is defined, as that will be
run by ``nose``, and a ``'__main__'`` block is defined too,
to allow running the script from the command line.

It is slightly non-standard, as :func:`NetlistTest.setUp()` and
:func:`NetlistTest.tearDown()` are called inside ``test()``, but this
was found to be an acceptable compromise between complexity and following
standard practices.

The script is meant to be run from the command line in case a regression
is detected by ``nose``, possibly with the aid of a debugger.
As such, the :func:`NetlistTest.tearDown()` function is not executed
in the ``'__main__'`` block, so that the test outputs are preserved for
inspection.

That said, the example file should be easy to understand and in most cases
a simple:

::

    :%s/mytest/<test_id>/g

in VIM - will suffice to generate your own script file. Just remember to save
to ``test_<test_id>.py``.

Overview of a typical test based on :class:`APITest`
\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"

Required files
^^^^^^^^^^^^^^

The main directory must contain:

- ``test_<test_id>.py``, the script executing the test,

- the reference data files for checking the pass/fail status of the test.
  These can be automatically generated, as it will be shown below.

Script file
'''''''''''

Again, it is probably easier to introduce the API test scripts with an example.

Below is a typical test script file:

::

    import ahkab
    from ahkab import ahkab, circuit, printing, devices, testing

    cli = False

    def test():
        \"\"\"Test docstring to be printed out by nose\"\"\"

        mycircuit = circuit.Circuit(title="Butterworth Example circuit", filename=None)

        ## define nodes
        gnd = mycircuit.get_ground_node()
        n1 = mycircuit.create_node('n1')
        n2 = mycircuit.create_node('n2')
        # ...

        ## add elements
        mycircuit.add_resistor(name="R1", n1="n1", n2="n2", value=600)
        mycircuit.add_inductor(name="L1", n1="n2", n2=gnd, value=15.24e-3)
        mycircuit.add_vsource("V1", n1="n1", n2=gnd, dc_value=5, ac_value=.5)
        # ...

        if cli:
            print(mycircuit)

        ## define analyses
        op_analysis = ahkab.new_op(outfile='<test_id>')
        ac_analysis = ahkab.new_ac(start=1e3, stop=1e5, points=100, outfile='<test_id>')
        # ...

        ## create a testbench
        testbench = testing.APITest('<test_id>', mycircuit,
                                    [op_analysis, ac_analysis],
                                    skip_on_travis=True, skip_on_pypy=True)

        ## setup and test
        testbench.setUp()
        testbench.test()

        ## this section is recommended. If something goes wrong, you may call the
        ## test from the cli and the plots to video in the following will allow
        ## for quick inspection
        if cli:
            ## re-run the test to grab the results
            r = ahkab.run(mycircuit, an_list=[op_analysis, ac_analysis])
            ## plot and save interesting data
            fig = plt.figure()
            plt.title(mycircuit.title + " - TRAN Simulation")
            plt.plot(r['tran']['T'], r['tran']['VN1'], label="Input voltage")
            plt.hold(True)
            plt.plot(r['tran']['T'], r['tran']['VN4'], label="output voltage")
            plt.legend()
            plt.hold(False)
            plt.grid(True)
            plt.ylabel('Step response')
            plt.xlabel('Time [s]')
            fig.savefig('tran_plot.png')
        else:
            ## don't forget to tearDown the testbench when under nose!
            testbench.tearDown()

    if __name__ == '__main__':
        import pylab as plt
        cli = True
        test()
        plt.show()

Once again, a function ``test()`` is defined, as that will be the
entry point of ``nose``, and a ``'__main__'`` block is defined as well,
to allow running the script from the command line.

Inside ``test()``, the circuit to be tested is defined, accessing the
``ahkab`` module directly, to set up elements, sources and analyses.
Directly calling :func:`ahkab.run()` is not necessary,
:func:`APITest.test()` will take care of that for you.

Notice how :func:`APITest.setUp()` and :func:`APITest.tearDown()` are
called inside ``test()``, as in the previous case.

The script is meant to be run from the command line in case a regression
is detected by ``nose``, possibly with the aid of a debugger.
As such, the :func:`APITest.tearDown()` function is not executed
in the ``'__main__'`` block, so that the test outputs are preserved for
inspection.

Additionally, plotting is performed if the test is directly run from
the command line.

In case non-standard simulation options are necessary, they can be set
as in the previous example.

Module reference
\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"

"""

from __future__ import (unicode_literals, absolute_import,
                        division, print_function)

import time
import os
import sys
import pickle
import unittest
from warnings import warn

try:
    from configparser import ConfigParser, NoOptionError
except ImportError:
    from ConfigParser import ConfigParser, NoOptionError

import numpy as np
import sympy

from scipy.interpolate import InterpolatedUnivariateSpline
from sympy.parsing.sympy_parser import parse_expr

from nose.tools import ok_, nottest
from nose.plugins.skip import SkipTest

from . import csvlib
from . import options
from . import py3compat
from . import pz
from . import results
from .ahkab import main, run


class _MyDict(dict):
    pass

@nottest
class NetlistTest(unittest.TestCase):
    """A class to run a netlist file and check the results against
    a pre-computed reference.

    **Parameters:**

    test_id : string
        The test id. For a netlist named ``"rc_network.ckt"``, this is
        to be set to ``"rc_network"``.

    er : float, optional
        Allowed relative error (applies to numeric results only).

    er : float, optional
        Allowed absolute error (applies to numeric results only).

    sim_opts : dict, optional
        A dictionary containing the options to be used for the test.

    verbose : int
        The verbosity level to be used in the test. From 0 (silent) to
        6 (verbose). Notice higher verbosity values usually result in
        higher coverage. Defaults to 6.
    """

    def __init__(self, test_id, er=1e-6, ea=1e-9, sim_opts=None, verbose=6):
        unittest.TestCase.__init__(self, methodName='test')
        self.test_id = test_id
        self.er = er
        self.ea = ea
        self.test.__func__.__doc__ = "%s simulation" % (test_id, )
        self.ref_data = {} # the reference results will be loaded here
        self._sim_opts = sim_opts if sim_opts is not None else {}
        self._reset_opts = {}
        self.verbose=verbose

    def _set_sim_opts(self, sim_opts):
        for opt in sim_opts.keys():
            if hasattr(options, opt):
                self._reset_opts.update({opt:getattr(options, opt)})
                setattr(options, opt, sim_opts[opt])
            else:
                raise ValueError("Option %s is not a valid option." % opt)

    def _reset_sim_opts(self):
        for opt in self._reset_opts:
            setattr(options, opt, self._reset_opts[opt])

    def setUp(self):
        """Set up the testbench."""
        # find the needed files wrt the WD
        # we may be called from <checkout-dir>/tests/<mytest>
        # or from <checkout-dir>/tests/
        # or from <checkout-dir>/
        wd = os.getcwd()
        if os.path.split(wd)[1] == self.test_id:
            self.reference_path = "."
        elif os.path.split(wd)[1] == 'tests':
            self.reference_path = os.path.join(wd, self.test_id)
        else:
            self.reference_path = os.path.join(wd, 'tests', self.test_id)

        if not os.path.isfile(os.path.join(self.reference_path,
                                           '%s.ini' % self.test_id)):
            raise IOError("Config file %s not found." %
                           os.path.join(self.reference_path,
                                        '%s.ini' % self.test_id))
        # read the test config from <test_id>.ini
        cp = ConfigParser()
        cp.read(os.path.join(self.reference_path, '%s.ini' % self.test_id))
        # skipping on TRAVIS-CI option for time-consuming tests
        self.skip = bool(int(cp.get('test', 'skip-on-travis')))
        if 'TRAVIS' in os.environ and self.skip:
            # skip even loading the references
            return
        # skipping on PYPY option for numeric tests
        # Do we have the optional skip-on-pypy entry?
        try:
            self.skip_on_pypy = bool(int(cp.get('test', 'skip-on-pypy')))
        except NoOptionError:
            # Default to skipping on PYPY
            self.skip_on_pypy = True
        if py3compat.PYPY and self.skip_on_pypy:
            # once again, skip even loading the references
            return

        assert self.test_id == cp.get('test', 'name')

        netlist = cp.get('test', 'netlist')
        self.netlist = os.path.join(self.reference_path, netlist)
        del netlist

        types = cp.get('test', 'type')
        self.types = [t.strip().replace(',', '').lower()
                      for t in types.split(',')]
        del types

        # reference files holding the reference results
        self.refs = {}
        for t in self.types:
            self.refs.update({t: os.path.join(self.reference_path,
                                 cp.get('test', t + '_ref'))})

        # files to be removed after the test is completed successfully
        self.rmfiles = []
        for i in self.types:
            if i == 'op':
                self.rmfiles.append(os.path.join(self.reference_path,
                                                 '%s.opinfo' %
                                                 self.test_id))
            self.rmfiles.append(os.path.join(self.reference_path,
                                                 '%s.%s' %
                                                 (self.test_id, i)))

        # Are we in a reference run?
        self.ref_run = False
        for i in list(self.refs.values()):
            self.ref_run = not os.path.isfile(i)
            if self.ref_run:
                print("RUNNING REFERENCE RUN - INVALID TEST!")
                break
        if not self.ref_run:
            self._load_references()

    def _load_references(self):
        for t, file_ref in list(self.refs.items()):
            if 'pickle' in file_ref:
                with open(file_ref, 'rb') as fp:
                    self.ref_data.update({t: pickle.load(fp)})
            else:
                data, headers, _, _ = csvlib.load_csv(file_ref, [], None, 0, verbose=0)
                res = _MyDict()
                if os.path.splitext(file_ref)[1][1:].lower() == 'ac':
                    res.update({headers[0]:data[0, :]})
                    for i, h in enumerate(headers):
                        if h[0] == h[-1] == '|':
                            pi = headers.index('arg('+h[1:-1]+')')
                            res.update({h[1:-1]:data[i, :]*np.exp(1j*data[pi, :])})
                        else:
                            continue
                else:
                    for i, h in enumerate(headers):
                        res.update({h: data[i, :]})
                res.x = headers[0]
                self.ref_data.update({t: res})

    def _run_test(self):
        # check whether we are on travis or not and skip if needed.
        # check whether we are running PYPY or not and skip if needed.
        if ('TRAVIS' in os.environ and self.skip) or (py3compat.PYPY and
            self.skip_on_pypy):
            self._reset_sim_opts()
            raise SkipTest
        # no reference runs with nose
        if sys.argv[0].endswith('nosetests') and self.ref_run:
            self._reset_sim_opts()
            raise SkipTest
        self._set_sim_opts(self._sim_opts)
        print("Running test... ", end="")
        start = time.time()
        res = main(filename=self.netlist,
                   outfile=os.path.join(self.reference_path, self.test_id),
                   verbose=self.verbose)
        stop = time.time()
        times = stop - start
        print("done.\nThe test took %f s" % times)
        return res

    def _check(self, res, ref):
        if hasattr(res, 'get_x'):
            x = res.get_x()
            for k in list(res.keys()):
                if np.all(res[k] == x):
                    continue
                elif np.any(np.iscomplex(res[k])) or np.any(np.iscomplex(ref[k])):
                    # Interpolate Re and Im of the results to compare.
                    x = x.reshape((-1, ))
                    refx = ref[ref.x].reshape((-1, ))
                    d1 = InterpolatedUnivariateSpline(x, np.real(res[k]).reshape((-1, )))
                    d2 = InterpolatedUnivariateSpline(refx, np.real(ref[k]).reshape((-1, )))
                    ok(d1(x), d2(x), rtol=self.er, atol=self.ea, msg=("Test %s FAILED (Re)" % self.test_id))
                    d1 = InterpolatedUnivariateSpline(x, np.imag(res[k]).reshape((-1, )))
                    d2 = InterpolatedUnivariateSpline(refx, np.imag(ref[k]).reshape((-1, )))
                    ok(d1(x), d2(x), rtol=self.er, atol=self.ea, msg=("Test %s FAILED (Im)" % self.test_id))
                else:
                    # Interpolate the results to compare.
                    x = x.reshape((-1, ))
                    refx = ref[ref.x].reshape((-1, ))
                    d1 = InterpolatedUnivariateSpline(x, np.real_if_close(res[k]).reshape((-1, )))
                    d2 = InterpolatedUnivariateSpline(refx, np.real_if_close(ref[k]).reshape((-1, )))
                    ok(d1(x), d2(x), rtol=self.er, atol=self.ea, msg=("Test %s FAILED" % self.test_id))
        elif isinstance(res, results.op_solution):
            for k in list(res.keys()):
                assert k in ref
                ok(res[k], ref[k], rtol=self.er, atol=self.ea, msg=("Test %s FAILED" % self.test_id))
        elif isinstance(res, results.pz_solution):
            # recover the reference signularities from Re/Im data
            ref_sing_keys = list(ref.keys())[:]
            ref_sing_keys.sort()
            assert len(ref_sing_keys) % 2 == 0
            ref_sing = [ref[ref_sing_keys[int(len(ref_sing_keys)/2) + k]] + ref[ref_sing_keys[k]]*1j \
                        for k in range(int(len(ref_sing_keys)/2))]
            ref_poles_num = len([k for k in ref.keys() if k[:4] == 'Re(p'])
            poles_ref, zeros_ref = ref_sing[:ref_poles_num], ref_sing[ref_poles_num:]
            assert len(poles_ref) == len(res.poles)
            pz._check_singularities(res.poles, poles_ref)
            assert len(zeros_ref) == len(res.zeros)
            pz._check_singularities(res.zeros, zeros_ref)
        else:
            if isinstance(res, list) or isinstance(res, tuple):
                for i, j in zip(res, ref):
                    self._check(i, j)
            elif res is not None:
                for k in list(res.keys()):
                    assert k in ref
                    if isinstance(res[k], dict): # hence ref[k] will be a dict too
                        self._check(res[k], ref[k])
                    elif isinstance(ref[k], sympy.Basic) and isinstance(res[k], sympy.Basic):
                        # get rid of assumptions. Evaluate only expression
                        rf = parse_expr(str(ref[k]))
                        rs = parse_expr(str(res[k]))
                        assert (rs == rf) or (sympy.simplify(rf/rs) == 1)
                    else:
                        assert res[k] == ref[k]

    def test(self):
        """Run the test."""
        res = self._run_test()
        if not self.ref_run:
            for t in list(res.keys()):
                ok_(t in self.ref_data, 'simulation %s not in the reference data' % t)
                print("Checking results for %s analysis..." % t)
                self._check(res[t], self.ref_data[t])
        else:
            for t, ref_file in list(self.refs.items()):
                if '.pickle' in ref_file:
                    with open(ref_file, 'wb') as fp:
                        pickle.dump(res[t], fp, protocol=2)
                else:
                    res_file = os.path.join(self.reference_path,
                                            '%s.%s' % (self.test_id, t))
                    os.rename(res_file, ref_file)

    def tearDown(self):
        """Remove temporary files - if needed."""
        if self.ref_run:
            pass
        else:
            for f in self.rmfiles:
                os.remove(f)
        self._reset_sim_opts()


@nottest
class APITest(unittest.TestCase):
    """A class to run a supplied circuit and check the results against
    a pre-computed reference.

    **Parameters:**

    test_id : string
        The test id.

    circ : circuit instance
        The circuit to be tested

    an_list : list of dicts
        A list of the analyses to be performed
    er : float, optional
        Allowed relative error (applies to numeric results only).

    er : float, optional
        Allowed absolute error (applies to numeric results only).

    sim_opts : dict, optional
        A dictionary containing the options to be used for the test.

    skip_on_travis : bool, optional
        Should we skip the test on Travis? Set to ``True`` for long tests.
        Defaults to ``False``.

    skip_on_pypy : bool, optional
        Should we skip the test on PYPY? Set to ``True`` for tests requiring
        libraries not supported by PYPY (eg. ``scipy``, ``matplotlib``).
        Defaults to ``True``, as most numeric tests will fail.
    """

    def __init__(self, test_id, circ, an_list, er=1e-6, ea=1e-9, sim_opts=None,
                 skip_on_travis=False, skip_on_pypy=True):
        unittest.TestCase.__init__(self, methodName='test')
        self.test_id = test_id
        self.er = er
        self.ea = ea
        self.test.__func__.__doc__ = "%s simulation" % (test_id, )
        self.ref_data = {} # the reference results will be loaded here
        self.skip = skip_on_travis
        self.skip_on_pypy = skip_on_pypy
        self.circ = circ
        self.an_list = an_list
        self._sim_opts = sim_opts if sim_opts is not None else {}
        self._reset_opts = {}
        self._set_sim_opts(self._sim_opts)
        self.res = None
        for an in an_list:
            if 'outfile' in an and self.test_id not in an['outfile']:
                warn("W: Analysis %s has outfile set to %s" %
                      (an['type'], an['outfile']))

    def _set_sim_opts(self, sim_opts):
        for opt in sim_opts.keys():
            if hasattr(options, opt):
                self._reset_opts.update({opt:getattr(options, opt)})
                setattr(options, opt, sim_opts[opt])
            else:
                raise ValueError("Option %s is not a valid option." % opt)

    def _reset_sim_opts(self):
        for opt in self._reset_opts:
            setattr(options, opt, self._reset_opts[opt])

    def setUp(self):
        """Set up the testbench"""
        # find the needed files wrt the WD
        # we may be called from <checkout-dir>/tests/<mytest>
        # or from <checkout-dir>/tests/
        # or from <checkout-dir>/
        wd = os.getcwd()
        if os.path.split(wd)[1] == self.test_id:
            self.reference_path = "."
        elif os.path.split(wd)[1] == 'tests':
            self.reference_path = os.path.join(wd, self.test_id)
        else:
            self.reference_path = os.path.join(wd, 'tests', self.test_id)

        if ('TRAVIS' in os.environ and self.skip) or (py3compat.PYPY and
            self.skip_on_pypy):
            # skip even loading the references
            return

        self.types = [a['type'] for a in self.an_list]

        # reference files holding the reference results
        self.refs = {}
        for t in self.types:
            self.refs.update({t: os.path.join(self.reference_path,
                                              self.test_id + '-ref' + '.'+ t)})

        # update the an_list with abs paths
        for i in range(len(self.an_list)):
            if 'outfile' in self.an_list[i] and \
               self.an_list[i]['outfile'] is not None and \
               not self.an_list[i]['outfile'] == 'stdout' and \
               not (len(self.an_list[i]['outfile']) > 5 and \
                    self.an_list[i]['outfile'][:4] == '/tmp/'):
                if not os.path.isabs(self.an_list[i]['outfile']):
                    self.an_list[i]['outfile'] = os.path.join(self.reference_path,
                                                       self.an_list[i]['outfile'])

        # files to be removed after the test is completed successfully
        self.rmfiles = []
        for an in self.an_list:
            if 'outfile' in an and \
               an['outfile'] is not None and \
               not an['outfile'] == 'stdout' and \
               not (len(an['outfile']) > 5 and an['outfile'][:4] == '/tmp/'):
                self.rmfiles.append(an['outfile'])
                if an['type'] == 'op':
                    self.rmfiles.append(an['outfile'] + 'info')

        # Are we in a reference run?
        self.ref_run = False
        for i in list(self.refs.values()):
            self.ref_run = not os.path.isfile(i)
            if self.ref_run:
                print("RUNNING REFERENCE RUN - INVALID TEST!")
                break
        if not self.ref_run:
            self._load_references()

    def _load_references(self):
        for t, file_ref in list(self.refs.items()):
            if '.symbolic' in file_ref:
                with open(file_ref, 'rb') as fp:
                    self.ref_data.update({t: pickle.load(fp)})
            else:
                data, headers, _, _ = csvlib.load_csv(file_ref, [], None, 0, verbose=0)
                res = _MyDict()
                if os.path.splitext(file_ref)[1][1:].lower() == 'ac':
                    res.update({headers[0]:data[0, :]})
                    for i, h in enumerate(headers):
                        if h[0] == h[-1] == '|':
                            pi = headers.index('arg('+h[1:-1]+')')
                            res.update({h[1:-1]:data[i, :]*np.exp(1j*data[pi, :])})
                        else:
                            continue
                else:
                    for i, h in enumerate(headers):
                        res.update({h: data[i, :]})
                res.x = headers[0] if not t == 'op' else None
                self.ref_data.update({t: res})

    def _run_test(self):
        if ('TRAVIS' in os.environ and self.skip) or (py3compat.PYPY and
            self.skip_on_pypy):
            self._reset_sim_opts()
            raise SkipTest
        print("Running test... ", end=' ')
        start = time.time()
        res = run(self.circ, self.an_list)
        stop = time.time()
        times = stop - start
        print("done.\nThe test took %f s" % times)
        return res

    def _check(self, res, ref):
        if hasattr(res, 'get_x'):
            x = res.get_x()
            for k in list(res.keys()):
                if np.all(res[k] == x):
                    continue
                elif np.any(np.iscomplex(res[k])) or np.any(np.iscomplex(ref[k])):
                    # Interpolate Re and Im of the results to compare.
                    x = x.reshape((-1, ))
                    refx = ref[ref.x].reshape((-1, ))
                    d1 = InterpolatedUnivariateSpline(x, np.real(res[k]).reshape((-1, )))
                    d2 = InterpolatedUnivariateSpline(refx, np.real(ref[k]).reshape((-1, )))
                    ok(d1(x), d2(x), rtol=self.er, atol=self.ea, msg=("Test %s FAILED (Re)" % self.test_id))
                    d1 = InterpolatedUnivariateSpline(x, np.imag(res[k]).reshape((-1, )))
                    d2 = InterpolatedUnivariateSpline(refx, np.imag(ref[k]).reshape((-1, )))
                    ok(d1(x), d2(x), rtol=self.er, atol=self.ea, msg=("Test %s FAILED (Im)" % self.test_id))
                else:
                    # Interpolate the results to compare.
                    x = x.reshape((-1, ))
                    refx = ref[ref.x].reshape((-1, ))
                    d1 = InterpolatedUnivariateSpline(x, np.real_if_close(res[k]).reshape((-1, )))
                    d2 = InterpolatedUnivariateSpline(refx, np.real_if_close(ref[k]).reshape((-1, )))
                    ok(d1(x), d2(x), rtol=self.er, atol=self.ea, msg=("Test %s FAILED" % self.test_id))
        elif isinstance(res, results.op_solution):
            for k in list(res.keys()):
                assert k in ref
                ok(res[k], ref[k], rtol=self.er, atol=self.ea, msg=("Test %s FAILED" % self.test_id))
        else:
            if isinstance(res, list) or isinstance(res, tuple):
                self._check(res[0], ref)
            elif res is not None:
                for k in list(res.keys()):
                    assert k in list(ref.keys())
                    if isinstance(res[k], dict): # hence ref[k] will be a dict too
                        self._check(res[k], ref[k])
                    elif isinstance(ref[k], sympy.Basic) and isinstance(res[k], sympy.Basic):
                        # get rid of assumptions. Evaluate only expression
                        rf = parse_expr(str(ref[k]))
                        rs = parse_expr(str(res[k]))
                        assert (rs == rf) or (sympy.simplify(rf/rs) == 1)
                    else:
                        assert res[k] == ref[k]

    def test(self):
        """Run the test."""
        res = self._run_test()
        if not self.ref_run:
            for t in list(res.keys()):
                ok_(t in self.ref_data, 'simulation %s not in the reference data')
                print("Checking results for %s analysis..." % t)
                self._check(res[t], self.ref_data[t])
        else:
            # move ref files into place
            for an in self.an_list:
                ref_file = self.refs[an['type']]
                if not os.path.isabs(an['outfile']):
                    res_file = os.path.join(self.reference_path, an['outfile'])
                else:
                    res_file = an['outfile']
                os.rename(res_file, ref_file)

    def tearDown(self):
        """Remove temporary files - if needed."""
        if self.ref_run:
            pass
        else:
            for f in self.rmfiles:
                os.remove(f)
        self._reset_sim_opts()

def ok(x, ref, rtol, atol, msg):
    try:
        assert np.allclose(x, ref, rtol=rtol, atol=atol)
    except AssertionError:
        print("REL: %g (max %g), ABS: %g (max %g)" % (max(abs(2*(x-ref)/(x+ref))), rtol, max(abs(x-ref)), atol))
        raise AssertionError(msg)

