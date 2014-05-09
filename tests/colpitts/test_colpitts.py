import time
import os, os.path
import pickle

import numpy
from scipy.interpolate import InterpolatedUnivariateSpline

from nose.tools import ok_, nottest, with_setup

import ahkab
from ahkab.ahkab import main

wd = os.getcwd()
if os.path.split(wd)[1] == 'ahkab':
	reference_path = os.path.join(wd, 'tests/colpitts')
elif os.path.split(wd)[1] == 'tests':
	reference_path = os.path.join(wd, 'colpitts')
else:
	reference_path = wd

ahkab_path = "ahkab"
er = 1e-6
ea = 1e-9


def _run_test(ref_run=False):
	netlist = os.path.join(reference_path, "colpitts.ckt")
	if not ref_run:
		data_file = os.path.join(reference_path, "colpitts")
	else:
		data_file = os.path.join(reference_path, "colpitts-ref")
	print "Running test... "
	start = time.time()
	main(filename=netlist, outfile=data_file, verbose=0)
	stop = time.time()
	times = stop - start
	print "Done. The test took %f s" % times
	data = numpy.loadtxt(data_file+".tran")
	return data, times

def teardown_func():
	for f in (os.path.join(reference_path, "colpitts.tran"), 
                  os.path.join(reference_path, "colpitts.opinfo")):
		os.remove(f)

@with_setup(None, teardown_func)
def test():
	"""Colpitts simulation"""
	ref_file = os.path.join(reference_path, 'colpitts-ref.tran')
	ref_run = not os.path.isfile(ref_file)

	if not ref_run:
		res = numpy.loadtxt(ref_file)
	else:
		print "RUNNING REFERENCE RUN - INVALID TEST!"

	res_new, time_new = _run_test(ref_run)
	
	# Interpolate the results to compare.
	d1 = InterpolatedUnivariateSpline(res_new[:, 0], res_new[:, 2])
	d2 = InterpolatedUnivariateSpline(res[:, 0], res[:, 2])

	ok_(numpy.allclose(d1(res[:, 0]), d1(res[:, 0]), rtol=er, atol=ea), "Test colpitts FAILED")

if __name__ == '__main__':
	test()
	ref_file = os.path.join(reference_path, 'colpitts-ref.tran')
	data = numpy.loadtxt(ref_file)
