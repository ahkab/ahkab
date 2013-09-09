import shooting
import bfpss
import options

def pss_analysis(**args):
	m = args.pop('method').lower()
	if m == options.SHOOTINGPSS:
		r = shooting.shooting(**args)
	elif m == options.BFPSS:
		r = bfpss.bfpss(**args)
	else:
		raise Exception, "Unknown PSS method %s" % m
	return r
