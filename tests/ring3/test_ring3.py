import time
import os, os.path

import numpy
from scipy.interpolate import InterpolatedUnivariateSpline

from nose.tools import ok_, nottest, with_setup

import ahkab
from ahkab.ahkab import main

ahkab_path = "ahkab"
er = 1e-6
ea = 1e-9

wd = os.getcwd()
if os.path.split(wd)[1] == 'ahkab':
	reference_path = os.path.join(wd, 'tests/ring3')
elif os.path.split(wd)[1] == 'tests':
	reference_path = os.path.join(wd, 'ring3')
else:
	reference_path = wd


def _run_test(ref_run=False):
	netlist = os.path.join(reference_path, "ring3.ckt")
	if not ref_run:
		data_file = os.path.join(reference_path, "ring3")
	else:
		data_file = os.path.join(reference_path, "ring3-ref")
	print "Running test... "
	start = time.time()
	main(filename=netlist, outfile=data_file, verbose=0)
	stop = time.time()
	times = stop - start
	print "Done. The test took %.3f s" % times
	data = numpy.loadtxt(data_file+".tran")
	return data, times

def teardown_func():
	for f in (os.path.join(reference_path, "ring3.tran"), 
	          os.path.join(reference_path, "ring3.opinfo")):
		os.remove(f)

@with_setup(None, teardown_func)
def test():
	"""Ring oscillator"""
	ref_file = os.path.join(reference_path, 'ring3-ref.tran')
	ref_run = not os.path.isfile(ref_file)

	if not ref_run:
		res = numpy.loadtxt(ref_file)
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
