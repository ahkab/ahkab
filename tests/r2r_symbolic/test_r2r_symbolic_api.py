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
import base64
import numpy

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
    mac = hex(uuid.getnode()).encode('ascii')
    #return base64.b64encode(hashlib.sha512(mac).hexdigest())
    hashlib.sha512(mac)
    hashlib.sha512(mac).hexdigest()
    hashlib.sha512(mac).hexdigest()#.decode('ascii')
    return base64.b64encode(hashlib.sha512(mac).hexdigest().encode('ascii'))

def load_boxid(filename):
    fp = open(filename, 'r')
    boxid = fp.read()
    fp.close()
    return boxid

def save_boxid(filename):
    fp = open(filename, 'w')
    boxid = get_boxid()
    fp.write(boxid.decode('ascii'))
    fp.flush()
    fp.close()

def check_boxid(filename):
    return get_boxid() == load_boxid(filename)

def _run_test(ref_run=False):
    MINNODES = 6
    MAXNODES = 14
    STEP = 1
    times = []

    x = list(range(max((2, MINNODES)), MAXNODES, STEP))
    for circuit_nodes in x:
        # build the circuit
        mycir = ahkab.Circuit('R2R symbolic test with %d nodes' %
                              circuit_nodes)
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
        out_th = 1./(2**(circuit_nodes - 1))
        assert .5*abs(out_th - out_test)/(out_th + out_test) < 1e-3

    x = numpy.array(x, dtype=numpy.int64)
    times = numpy.array(times, dtype=numpy.float64)
    if ref_run:
        numpy.savetxt(os.path.join(reference_path, 'r2r_symbolic_ref.csv'),
                      numpy.concatenate((x.reshape((-1, 1)),
                                         times.reshape((-1, 1))),
                                        axis=1), delimiter=',')
        save_boxid(os.path.join(reference_path, 'r2r_symbolic_ref.boxid'))
    return x, times


def test():
    """R-2R ladder symbolic speed test"""

    # we do not want to execute this on Travis.
    #if 'TRAVIS' in os.environ:
    #    # we skip the test. Travis builders are awfully slow
    #    raise SkipTest

    csv_file = os.path.join(reference_path, 'r2r_symbolic_ref.csv')
    boxid_file = os.path.join(reference_path, 'r2r_symbolic_ref.boxid')
    ref_run = not (os.path.isfile(csv_file) and check_boxid(boxid_file))
    if not ref_run:
        print("Running test...")
        d = numpy.loadtxt(csv_file, delimiter=',')
        x, times = d[:, 0].reshape((-1,)), d[:, 1].reshape((-1,))
        x = numpy.array(x, dtype=numpy.int64)
        x_new, times_new = _run_test(ref_run)
    elif ref_run:
        if sys.argv[0].endswith('nosetests'):
            raise SkipTest
        print("RUNNING REFERENCE TEST - RESULTS INVALID!")
        x, times = _run_test(ref_run)
        print("RUNNING REFERENCE TEST - RESULTS INVALID!")

    if stand_alone_exec:
        image_file = os.path.join(reference_path, "r2r_symbolic.txt")
        image_file_obj = open(image_file, 'w')
        if not ref_run:
            plot_comparison(x, times, times_new, label1='Measured - REF',
                            label2='Measured - NEW', fileobj=image_file_obj)
        else:
            plot_comparison(x, times, [], label1='Measured - REF',
                            label2='Measured - NEW', fileobj=image_file_obj)
        image_file_obj.flush()
        image_file_obj.close()

    # last thing we do is we check for timing issues.
    if not ref_run:
        assert numpy.max(times_new - times) < REGRESSION_TIME
        assert numpy.sum(times_new) > 3  # if we're that fast, something's off

def plot_comparison(x, y1, y2, label1=None, label2=None, fileobj=None):
    """Plot by printing histograms to screen

    **Parameters:**

    x : array-like
        The abscissa values for ``y1`` and ``y2``.
    y1 : array-like
        The ordinata values for the first histogram.
    y2 : array-like
        The ordinata values for the second data set. It may be set to a
        zero-length array if only one data set is to be plotted.
    label1 : str, optional
        The label for the first data set.
    label2 : str, optional
        The label for the second data set.
    fileobj : file object, optional
        If you wish to redirect the output to a special file, you may.
        Otherwise, ``sys.stdout`` will be employed.

    """
    if fileobj is None:
        fileobj = sys.stdout
    maxy2 = max(y2) if len(y2) else 1.
    maxy1 = max(y1) if len(y1) else 1.
    maxv = float(max([1., maxy1, maxy2]))
    scale = lambda x: int(x/maxv*80)
    unscale = lambda x: x*maxv/80
    # top x
    print('         ', end=' ', file=fileobj)
    for i in range(9):
        print('%5.3f    ' % unscale(i*10), end=' ', file=fileobj)
    print(file=fileobj)
    print('          |        ', end=' ', file=fileobj)
    for i in range(1, 9):
        print('|        ', end=' ', file=fileobj)
    print(file=fileobj)
    print('          +'+'-'*81+'>', file=fileobj)
    # data
    print(('          |'), file=fileobj)
    if len(y2):
        for xi, y1i, y2i in zip(x, y1, y2):
            print(('% 9.3f +' % xi) + '-'*(scale(y1i)-1)+'o', file=fileobj)
            print(('          |') + '='*(scale(y2i)-1)+'*', file=fileobj)
            print(('          |'), file=fileobj)
    else:
        for xi, y1i in zip(x, y1):
            print(('% 9.3f +' % xi) + '-'*(scale(y1i)-1)+'o', file=fileobj)
            print(('          |'), file=fileobj)
    # bottom x
    print('          +'+'-'*81+'>', file=fileobj)
    print('          |        ', end=' ', file=fileobj)
    for i in range(1, 9):
        print('|        ', end=' ', file=fileobj)
    print(file=fileobj)
    print('         ', end=' ', file=fileobj)
    for i in range(9):
        print('%5.3f    ' % unscale(i*10), end=' ', file=fileobj)
    if label1 or label2:
        print('\n\n   LEGEND:', file=fileobj)
    if label1:
        print('   --o %s' % label1, file=fileobj)
    if label2 and len(y2):
        print('   ==* %s' % label2, file=fileobj)

if __name__ == '__main__':
    stand_alone_exec = True
    test()
