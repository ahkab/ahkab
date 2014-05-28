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
import io
import sys
import hashlib
import uuid
import scipy.io
import subprocess
import base64
import scipy, scipy.optimize
import numpy
import matplotlib
matplotlib.use('Agg')
import pylab

from nose.plugins.skip import SkipTest

stand_alone_exec = False
REGRESSION_TIME = 3 #s - if the test takes 3s additional time it's considered
                    # a regression and the test fails

# find out the wd so we can find the files loc.
wd = os.getcwd()
if os.path.split(wd)[1] == 'ahkab':
    reference_path = os.path.join(wd, 'tests/r2r')
elif os.path.split(wd)[1] == 'tests':
    reference_path = os.path.join(wd, 'r2r')
else:
    reference_path = wd

# checking for time regression on different machines makes little sense
# this is an easy -- and imperfect -- way to check whether ref and
# current execution happen on the same box
def get_boxid():
    mac = hex(uuid.getnode())
    return base64.b64encode(hashlib.sha512(mac).hexdigest())

def check_boxid(mat_file):
    return scipy.io.loadmat(mat_file)['boxid']

# fit the data
fitfunc = lambda p, x: p[0]*x**3 + p[1]*x**2 + p[2] # Target function
errfunc = lambda p, x, y: (fitfunc(p, x) - y) # Distance to the target function

def fit(y, x):
    fit = [3.1, 1.1, .2]
    fit, success = scipy.optimize.leastsq(errfunc, fit, args=(x, y), maxfev=100000)
    print(fit, success)
    return fit

def _run_test(ref_run=False, verbose=False):
    title = "* MATRIX SIZE TEST: R-2R ladder\n"
    start_block = "V1 1 0 type=vdc vdc=10e3\n"
    middle_block1 = "R%(block)dh %(input)d %(output)d 1.3k\n"
    middle_block2 = "R%(block)dv %(output)d 0 2.6k\n"
    last_block = "R%(block)dve %(input)d 0 2.6k\n"
    analysis = ".op\n"
    
    MAXNODES = 2000
    STEP = 50

    filename = "r2r_%(nodes)s.ckt"
    times = []

    for circuit_nodes in range(0, MAXNODES, STEP):
        tmp_filename = os.path.join(reference_path, filename % {'nodes':circuit_nodes})
        if circuit_nodes < 2: continue
        fp = open(tmp_filename, 'w')
        fp.write(title)
        fp.write(start_block)
        for n in range(1, circuit_nodes):
            fp.write(middle_block1 % {'block':n, 'input':n, 'output':n+1})
            fp.write(middle_block2 % {'block':n, 'input':n, 'output':n+1})
        fp.write(last_block % {'block':circuit_nodes, 'input':circuit_nodes})
        fp.write(analysis)
        fp.close()
        start = time.time()
        outstr = subprocess.check_output(["ahkab", "-v 0", tmp_filename])
        stop = time.time()
        times.append((stop - start))
        if verbose:
            with open(tmp_filename, 'r') as fp:
                print("".join(fp.readlines()))
            print(outstr)
        print("Solving with %d nodes took %f s" % (circuit_nodes, times[-1]))
        os.remove(tmp_filename)

    x = list(range(0, MAXNODES, STEP))[1:]
    x = numpy.array(x, dtype='int64')
    times = numpy.array(times, dtype='float64')
    if ref_run:
        scipy.io.savemat(os.path.join(reference_path, 'r2r.mat'), 
                         {'x':x, 'times':times, 'boxid':get_boxid()})
    return x, times

def test():
    """R-2R ladder speed test"""

    # we do not want to execute this on Travis.
    if 'TRAVIS' in os.environ:
        # we skip the test. Travis builders are awfully slow
        raise SkipTest 

    mat_file = os.path.join(reference_path, 'r2r.mat') 
    ref_run = not (os.path.isfile(mat_file) and check_boxid(mat_file))
    if not ref_run:
        print("Running test...")
        d = scipy.io.loadmat(mat_file)
        x, times = d['x'].reshape((-1,)), d['times'].reshape((-1,))
        x = numpy.array(x, dtype='int64')
        x_new, times_new = _run_test(ref_run)
        assert numpy.max(times_new - times) < REGRESSION_TIME
        assert numpy.sum(times_new) > 3 # if we're that fast, something's off
    elif ref_run:
        if sys.argv[0].endswith('nosetests'):
            raise SkipTest
        print("RUNNING REFERENCE TEST - RESULTS INVALID!")
        x, times = _run_test(ref_run)
        print("RUNNING REFERENCE TEST - RESULTS INVALID!")

    if stand_alone_exec:
        image_file = os.path.join(reference_path, "r2r.png")
        pylab.figure()
        pylab.hold(True)
        pylab.title("Total times vs number of equations")
        pylab.plot(x, times, 'ko', label='Measured - REF')
        p = fit(times, x)
        xf = numpy.arange(x[-1])
        pylab.plot(xf, fitfunc(p, xf), 'k-', label=('Fit: $y = %.3e\ x^3 + %.3e\ x^2 + %.1f$' % tuple(p.tolist())))
        if not ref_run:
            pylab.plot(x_new, times_new, 'go', label='Measured - NEW')
            p_new = fit(times_new, x_new)
            xf_new = numpy.arange(x_new[-1])
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
