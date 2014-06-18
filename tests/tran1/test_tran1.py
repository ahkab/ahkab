from __future__ import unicode_literals, print_function
from ahkab.testing import NetlistTest
from ahkab import options
options.plotting_show_plots = False

def myoptions(TM):
    sim_opts = dict(default_tran_method=TM)
    return sim_opts

def test():
    print("")
    tms = ('IMPLICIT_EULER', 'TRAP', 'GEAR3', 'GEAR4', 'GEAR5', 'GEAR6')
    for TMI, TM in enumerate(tms):
        print('%d/%d Testing %s...' % (TMI+1, len(tms), TM))
        nt = NetlistTest('tran1', sim_opts=myoptions(TM))
        nt.setUp()
        nt.test()
        nt.tearDown()

test.__doc__ = "Simple TRAN methods test"

if __name__ == '__main__':
    nt = NetlistTest('tran1')
    nt.setUp()
    nt.test()
