# -*- coding: utf-8 -*-
# test_r2r.py
# R2r ladder timed test
# Copyright 2012 Giuseppe Venturini
#
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
from __future__ import print_function, division, unicode_literals
import time
import os
import sys
import hashlib
import uuid
import scipy.io
import base64
import scipy
import scipy.optimize
import numpy
import matplotlib
matplotlib.use('Agg')
import pylab

from nose.plugins.skip import SkipTest

import ahkab
from ahkab import py3compat

stand_alone_exec = False
# if the test takes 3s additional time it's considered
# a regression and the test fails
REGRESSION_TIME = 3  # seconds

# find out the wd so we can find the files loc.
wd = os.getcwd()
if os.path.split(wd)[1] == 'ahkab':
    reference_path = os.path.join(wd, 'tests/r2r_symbolic')
elif os.path.split(wd)[1] == 'tests':
    reference_path = os.path.join(wd, 'r2r_symbolic')
else:
    reference_path = wd

# checking for time regression on different machines makes little sense
# this is an easy -- and imperfect -- way to check whether ref and
# current execution happen on the same box


def get_boxid():
    mac = hex(uuid.getnode())
    return base64.b64encode(hashlib.sha512(mac).hexdigest())


def check_boxid(mat_file):
    return get_boxid() == scipy.io.loadmat(mat_file)['boxid']

# fit the data
fitfunc = lambda p, x: p[0] * x**3 + p[1] * x**2 + p[2]  # Target function
# Distance to the target function
errfunc = lambda p, x, y: (fitfunc(p, x) - y)


def fit(y, x):
    fit = [3.1, 1.1, .2]
    fit, success = scipy.optimize.leastsq(
        errfunc, fit, args=(x, y), maxfev=100000)
    print(fit, success)
    return fit


def _run_test(ref_run=False):
    MAXNODES = 14
    STEP = 1

    filename = "r2r_symbolic_%(nodes)s.ckt"
    times = []

    for circuit_nodes in range(2, MAXNODES, STEP):
        # build the circuit
        mycir = ahkab.Circuit(
            'R2R symbolic test with %d nodes' % circuit_nodes)
        n1 = '1'
        gnd = mycir.gnd
        mycir.add_vsource('VS', n1, gnd, dc_value=10e3)
        subs = {}
        for n in range(1, circuit_nodes):
            n1 = str(n)
            n2 = str(n + 1)
            mycir.add_resistor("R%dh1" % n, n1, n2, value=2.6e3)
            mycir.add_resistor("R%dh2" % n, n1, n2, value=2.6e3)
            mycir.add_resistor("R%dv" % n, n2, gnd, value=2.6e3)
            subs.update({"R%dh1" % n: 'R', "R%dh2" % n: 'R', "R%dv" % n: 'R'})
        n1 = str(circuit_nodes)
        mycir.add_resistor("R%dve" % circuit_nodes, n1, gnd, value=2.6e3)
        subs.update({"R%dve" % circuit_nodes: 'R'})
        # define analysis
        s = ahkab.new_symbolic(subs=subs)
        start = time.time()
        r = ahkab.run(mycir, s)['symbolic'][0]
        stop = time.time()
        times.append((stop - start))
        print("Solving with %d nodes took %f s" % (circuit_nodes, times[-1]))
        # check the values too
        VS = r.as_symbol('VS')
        out_test = r['V' + str(circuit_nodes)] / VS
        out_th = 1. / (2**(circuit_nodes - 1))
        assert .5 * abs(out_th - out_test) / (out_th + out_test) < 1e-3

    x = list(range(2, MAXNODES, STEP))
    x = numpy.array(x, dtype='int64')
    times = numpy.array(times, dtype='float64')
    if ref_run:
        scipy.io.savemat(os.path.join(reference_path, 'r2r_symbolic.mat'),
                         {'x': x, 'times': times, 'boxid': get_boxid()})
    return x, times


def test():
    """R-2R ladder symbolic speed test"""

    # we do not want to execute this on Travis.
    if 'TRAVIS' in os.environ or py3compat.PYPY: 
        # we skip the test. Travis builders are awfully slow
        # also no execution on PYPY for now, we need scipy for the
        # .mat-files interface 
        raise SkipTest

    mat_file = os.path.join(reference_path, 'r2r_symbolic.mat')
    ref_run = not (os.path.isfile(mat_file) and check_boxid(mat_file))
    if not ref_run:
        print("Running test...")
        d = scipy.io.loadmat(mat_file)
        x, times = d['x'].reshape((-1,)), d['times'].reshape((-1,))
        x = numpy.array(x, dtype='int64')
        x_new, times_new = _run_test(ref_run)
        assert numpy.max(times_new - times) < REGRESSION_TIME
        assert numpy.sum(times_new) > 3  # if we're that fast, something's off
    elif ref_run:
        if sys.argv[0].endswith('nosetests'):
            raise SkipTest
        print("RUNNING REFERENCE TEST - RESULTS INVALID!")
        x, times = _run_test(ref_run)
        print("RUNNING REFERENCE TEST - RESULTS INVALID!")

    if stand_alone_exec:
        image_file = os.path.join(reference_path, "r2r_symbolic.png")
        pylab.figure()
        pylab.hold(True)
        pylab.title("Total times vs number of equations")
        pylab.plot(x, times, 'ko', label='Measured - REF')
        p = fit(times, x)
        xf = numpy.arange(2, x[-1] + 1)
        pylab.plot(xf, fitfunc(p, xf), 'k-',
                   label=('Fit: $y = %.3e\ x^3 + %.3e\ x^2 + %.1f$' % tuple(p.tolist())))
        if not ref_run:
            pylab.plot(x_new, times_new, 'go', label='Measured - NEW')
            p_new = fit(times_new, x_new)
            xf_new = numpy.arange(2, x_new[-1] + 1)
            pylab.plot(xf_new, fitfunc(p_new, xf_new), 'k-',
                       label=('Fit: $y = %.3e\ x^3 + %.3e\ x^2 + %.1f$ - NEW' % tuple(p_new.tolist())))
        pylab.xlabel("Number of equations []")
        pylab.ylabel("Time to convergence [s]")
        pylab.grid(True)
        pylab.legend(loc=0)
        pylab.savefig(image_file, dpi=90, format='png')

if __name__ == '__main__':
    stand_alone_exec = True
    test()
