import shooting
import bfpss
import options

# .SHOOTING PERIOD=n [points=n step=n autonomous=bool]

specs = {'pss':{'tokens':({
                          'label':'period',
                          'pos':0,
                          'type':float,
                          'needed':True,
                          'dest':'period',
                          'default':None
                         },
                         {
                          'label':'points',
                          'pos':None,
                          'type':float,
                          'needed':True,
                          'dest':'points',
                          'default':None
                         },
                         {
                          'label':'step',
                          'pos':None,
                          'type':float,
                          'needed':False,
                          'dest':'step',
                          'default':None
                         },
                         {
                          'label':'autonomous',
                          'pos':None,
                          'type':bool,
                          'needed':False,
                          'dest':'uic',
                          'default':'0'
                         },
                         {
                          'label':'method',
                          'pos':None,
                          'type':str,
                          'needed':False,
                          'dest':'method',
                          'default':options.SHOOTINGPSS
                         }
                        )
               }
           }

def pss_analysis(*largs, **args):
	m = args.pop('method').lower()
	if m == options.SHOOTINGPSS:
		r = shooting.shooting(*largs, **args)
	elif m == options.BFPSS:
		r = bfpss.bfpss(*largs, **args)
	else:
		raise Exception, "Unknown PSS method %s" % m
	return r
