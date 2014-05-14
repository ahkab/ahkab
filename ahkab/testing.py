import time
import os
import pickle
import unittest

from ConfigParser import ConfigParser

import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline

from nose.tools import ok_, nottest

from . import csvlib
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
    """

    def __init__(self, test_id, er=1e-6, ea=1e-9):
        unittest.TestCase.__init__(self, methodName='test')
        self.test_id = test_id
        self.er = er
        self.ea = ea
        self.test.__func__.__doc__ = "%s simulation" % (test_id, )
        self.ref_data = {} # the reference results will be loaded here

    def setUp(self):
        # find the needed files wrt the WD
        # we may be called from ahkab/tests/<mytest>
        # or from tests/<mytest>
        # or from <mytest>
        wd = os.getcwd()
        if os.path.split(wd)[1] == 'ahkab':
            self.reference_path = os.path.join(wd, 'tests', self.test_id)
        elif os.path.split(wd)[1] == 'tests':
            self.reference_path = os.path.join(wd, self.test_id)
        else:
            self.reference_path = "."

        # read the test config from <test_id>.ini
        cp = ConfigParser()
        cp.read(os.path.join(self.reference_path, '%s.ini' % self.test_id))
        self.skip = bool(cp.get('test', 'skip-on-travis'))
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
        for i in self.refs.values():
            self.ref_run = not os.path.isfile(i)
            if self.ref_run:
                print "RUNNING REFERENCE RUN - INVALID TEST!"
                break
        if not self.ref_run:
            self.load_references()

    def load_references(self):
        for t, file_ref in self.refs.items():
            if 'pickle' in file_ref:
                with open(file_ref, 'r') as fp:
                    self.ref_data.update({t: pickle.load(fp)})
            else:
                data, headers, _, _ = csvlib.load_csv(file_ref, [], None, 0L, verbose=0)
                res = _MyDict()
                for i, h in enumerate(headers):
                    res.update({h: data[i, :]})
                res.x = headers[0]
                self.ref_data.update({t: res})

    def _run_test(self):
        print "Running test... ",
        start = time.time()
        res = main(filename=self.netlist,
                   outfile=os.path.join(self.reference_path, self.test_id),
                   verbose=0)
        stop = time.time()
        times = stop - start
        print "done.\nThe test took %f s" % times
        return res

    def _check(self, res, ref):
        if hasattr(res, 'get_x'):
            x = res.get_x()
            for k in res.keys():
                if res[k] is x:
                    continue
                else:
                    # Interpolate the results to compare.
                    d1 = InterpolatedUnivariateSpline(x.reshape((-1, )), res[k].reshape((-1, )))
                    d2 = InterpolatedUnivariateSpline(ref[ref.x].reshape((-1, )), res[k].reshape((-1, )))
                    ok_(np.allclose(d1(x.reshape((-1, ))), d2(x.reshape((-1, ))), rtol=self.er, atol=self.ea), "Test %s FAILED" % self.test_id)
        elif isinstance(res, results.op_solution):
            for k in res.keys():
                assert k in ref
                ok_(np.allclose(res[k], ref[k], rtol=self.er, atol=self.ea), "Test %s FAILED" % self.test_id)
        else:
            if isinstance(res, list) or isinstance(res, tuple):
                for i, j in zip(res, ref):
                    self._check(i, j)
            elif res is not None:
                for k in res.keys():
                    assert k in ref
                    assert res[k] == ref[k]

    def test(self):
        res = self._run_test()
        if not self.ref_run:
            ok_(set(list(res.keys())) == set(list(self.ref_data.keys())),
                "Reference and test data have a different number of nodes")
            for t in res.keys():
                ok_(t in self.ref_data, 'simulation %s not in the reference data')
                print "Checking results for %s analysis..." % t
                self._check(res[t], self.ref_data[t])
        else:
            for t, ref_file in self.refs.items():
                if '.pickle' in ref_file:
                    with open(ref_file, 'w') as fp:
                        pickle.dump(res[t], fp)
                else:
                    res_file = os.path.join(self.reference_path, 
                                            '%s.%s' % (self.test_id, t))
                    os.rename(res_file, ref_file)

    def tearDown(self):
        if self.ref_run:
            pass
        else:
            for f in self.rmfiles:
                os.remove(f)


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

    skip_on_travis : bool, optional
        Should we skip the test on Travis? Set to ``True`` for long tests
    """

    def __init__(self, test_id, circ, an_list, er=1e-6, ea=1e-9, skip_on_travis=False):
        unittest.TestCase.__init__(self, methodName='test')
        self.test_id = test_id
        self.er = er
        self.ea = ea
        self.test.__func__.__doc__ = "%s simulation" % (test_id, )
        self.ref_data = {} # the reference results will be loaded here
        self.skip = skip_on_travis
        self.circ = circ
        self.an_list = an_list

    def setUp(self):
        # find the needed files wrt the WD
        # we may be called from ahkab/tests/<mytest>
        # or from tests/<mytest>
        # or from <mytest>
        wd = os.getcwd()
        if os.path.split(wd)[1] == 'ahkab':
            self.reference_path = os.path.join(wd, 'tests', self.test_id)
        elif os.path.split(wd)[1] == 'tests':
            self.reference_path = os.path.join(wd, self.test_id)
        else:
            self.reference_path = "."

        self.types = [a['type'] for a in self.an_list]

        # reference files holding the reference results
        self.refs = {}
        for t in self.types:
            self.refs.update({t: os.path.join(self.reference_path,
                                              self.test_id + '_ref' + '.'+ t)})

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
        for i in self.refs.values():
            self.ref_run = not os.path.isfile(i)
            if self.ref_run:
                print "RUNNING REFERENCE RUN - INVALID TEST!"
                break
        if not self.ref_run:
            self.load_references()

    def load_references(self):
        for t, file_ref in self.refs.items():
            if '.symbolic' in file_ref:
                with open(file_ref, 'rb') as fp:
                    self.ref_data.update({t: pickle.load(fp)})
            else:
                data, headers, _, _ = csvlib.load_csv(file_ref, [], None, 0L, verbose=0)
                res = _MyDict()
                for i, h in enumerate(headers):
                    res.update({h: data[i, :]})
                res.x = headers[0] if not t == 'op' else None
                self.ref_data.update({t: res})

    def _run_test(self):
        print "Running test... ",
        start = time.time()
        res = run(self.circ, self.an_list)
        stop = time.time()
        times = stop - start
        print "done.\nThe test took %f s" % times
        return res

    def _check(self, res, ref):
        if hasattr(res, 'get_x'):
            x = res.get_x()
            for k in res.keys():
                if res[k] is x:
                    continue
                else:
                    # Interpolate the results to compare.
                    d1 = InterpolatedUnivariateSpline(x.reshape((-1, )), res[k].reshape((-1, )))
                    d2 = InterpolatedUnivariateSpline(ref[ref.x].reshape((-1, )), res[k].reshape((-1, )))
                    ok_(np.allclose(d1(x.reshape((-1, ))), d2(x.reshape((-1, ))), rtol=self.er, atol=self.ea), "Test %s FAILED" % self.test_id)
        elif isinstance(res, results.op_solution):
            for k in res.keys():
                assert k in ref
                ok_(np.allclose(res[k], ref[k], rtol=self.er, atol=self.ea), "Test %s FAILED" % self.test_id)
        else:
            if isinstance(res, list) or isinstance(res, tuple):
                for i, j in zip(res, ref):
                    self._check(i, j)
            elif res is not None:
                for k in res.keys():
                    assert k in ref
                    assert res[k] == ref[k]

    def test(self):
        res = self._run_test()
        if not self.ref_run:
            ok_(set(list(res.keys())) == set(list(self.ref_data.keys())),
                "Reference and test data have a different number of nodes")
            for t in res.keys():
                ok_(t in self.ref_data, 'simulation %s not in the reference data')
                print "Checking results for %s analysis..." % t
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
        if self.ref_run:
            pass
        else:
            for f in self.rmfiles:
                os.remove(f)
