import time
import os, os.path
import pickle
import subprocess
import scipy, scipy.optimize
import numpy
import matplotlib
matplotlib.use('Agg')
import pylab

stand_alone_exec = False

# fit the code
fitfunc = lambda p, x: p[0]*x**3 + p[1]*x**2 + p[2] # Target function
errfunc = lambda p, x, y: (fitfunc(p, x) - y) # Distance to the target function

def fit(y, x):
        fit = [3.1, 1.1, .2]
        fit, success = scipy.optimize.leastsq(errfunc, fit, args=(x, y), maxfev=100000)
        print fit, success
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
		if circuit_nodes < 2: continue
		fp = open(filename % {'nodes':circuit_nodes}, 'w')
		fp.write(title)
		fp.write(start_block)
		for n in range(1, circuit_nodes):
			fp.write(middle_block1 % {'block':n, 'input':n, 'output':n+1})
			fp.write(middle_block2 % {'block':n, 'input':n, 'output':n+1})
		fp.write(last_block % {'block':circuit_nodes, 'input':circuit_nodes})
		fp.write(analysis)
		fp.close()
		start = time.time()
		outstr = subprocess.check_output(["ahkab", "-v 0", filename % {'nodes':circuit_nodes}])
		stop = time.time()
		times.append((stop-start))
		if verbose:
			with open(filename % {'nodes':circuit_nodes}, 'r') as fp:
				print "".join(fp.readlines())
			print outstr
		print "Solving with %d nodes took %f s" % (circuit_nodes, times[-1])
		os.remove(filename % {'nodes':circuit_nodes})

	x = range(0, MAXNODES, STEP)[1:]
	x = numpy.array(x, dtype='int64')
	times = numpy.array(times, dtype='float64')
	if ref_run:
		with open('r2r.pickle', 'w') as fp:
			pickle.dump((x, times), fp)
	return x, times

def test():
	ref_run = not os.path.isfile('r2r.pickle')
	if not ref_run:
		x, times = pickle.load(open('r2r.pickle', 'r'))
		x = numpy.array(x, dtype='int64')
		x_new, times_new = _run_test(ref_run)
	else:
		print "RUNNING REFERENCE TEST - RESULTS INVALID!"
		x, times = _run_test(ref_run)
		print "RUNNING REFERENCE TEST - RESULTS INVALID!"

	if stand_alone_exec:
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
			pylab.plot(xf_new, fitfunc(p_new, xf_new), 'k-', label=('Fit: $y = %.3e\ x^3 + %.3e\ x^2 + %.1f$ - NEW' % tuple(p_new.tolist())))
		pylab.xlabel("Number of equations []")
		pylab.ylabel("Time to convergence [s]")
		pylab.grid(True)
		pylab.legend(loc=0)
		pylab.savefig("r2r.png", dpi=90, format='png')

if __name__ == '__main__':
	stand_alone_exec = True
	test()
