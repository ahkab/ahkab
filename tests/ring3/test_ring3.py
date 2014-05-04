import time
import os, os.path
import subprocess

import numpy
from scipy.interpolate import InterpolatedUnivariateSpline

from nose.tools import ok_, nottest, with_setup

ahkab_path = "ahkab"
er = 1e-6
ea = 1e-9


def _run_test(ref_run=False):
	netlist = "ring3.ckt"
	data_file = "ring3" if not ref_run else "ring3-ref"
	print "Running test... "
	start = time.time()
	proc = subprocess.Popen([ahkab_path, "-v", "0", "-o", data_file, netlist])
	proc.communicate()
	stop = time.time()
	times = stop-start
	print "Done. The test took %f s" % times
	data = numpy.loadtxt(data_file+".tran")
	return data, times

def teardown_func():
	for f in ("ring3.tran", "ring3.op"):
		os.remove(f)

@with_setup(None, teardown_func)
def test():
	ref_run = not os.path.isfile('ring3-ref.tran')

	if not ref_run:
		res = numpy.loadtxt("ring3-ref.tran")
	else:
		print "RUNNING REFERENCE RUN - INVALID TEST!"

	res_new, time_new = _run_test(ref_run)
	assert res_new.shape[1] == res.shape[1] # same number of unknowns

	# Interpolate the results to compare.
	d1s = []
	d2s = []
	for i in range(1, res_new.shape[1]):
		d1s += [InterpolatedUnivariateSpline(res_new[:, 0], res_new[:, i])]
        	d2s += [InterpolatedUnivariateSpline(res[:, 0], res[:, i])]

	for d1, d2 in zip(d1s, d2s): # we rely on the solutions ORDER. Not good, not good...
		ok_(numpy.allclose(d1(res[:, 0]), d2(res[:, 0]), rtol=er, atol=ea), "Test RING3 FAILED")

if __name__ == '__main__':
	test()
	teardown_func()
