import time
import os, os.path
import pickle
import subprocess
import numpy
from scipy.interpolate import InterpolatedUnivariateSpline
from nose.tools import ok_, nottest, with_setup

ahkab_path = "ahkab"
er = 1e-6
ea = 1e-9


def _run_test(ref_run=False):
	netlist = "colpitts.ckt"
	data_file = "colpitts" if not ref_run else "colpitts-ref"
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
	for f in ("colpitts.tran", "colpitts.op"):
		os.remove(f)

@with_setup(None, teardown_func)
def test():
	ref_run = not os.path.isfile('colpitts-ref.tran')

	if not ref_run:
		res = numpy.loadtxt("colpitts-ref.tran")
	else:
		print "RUNNING REFERENCE RUN - INVALID TEST!"

	res_new, time_new = _run_test(ref_run)
	
	# Interpolate the results to compare.
	d1 = InterpolatedUnivariateSpline(res_new[:, 0], res_new[:, 2])
	d2 = InterpolatedUnivariateSpline(res[:, 0], res[:, 2])

	ok_(numpy.allclose(d1(res[:, 0]), d1(res[:, 0]), rtol=er, atol=ea), "Test colpitts FAILED")

if __name__ == '__main__':
	test()
	data = numpy.loadtxt("colpitts-ref.tran")
