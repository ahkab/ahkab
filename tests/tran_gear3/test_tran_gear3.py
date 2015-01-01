from __future__ import unicode_literals, print_function
from ahkab.testing import NetlistTest
from ahkab import options
options.plotting_show_plots = False

def test():
    nt = NetlistTest('tran_gear3', ea=1e-3, er=1e-3)
    nt.setUp()
    nt.test()
    nt.tearDown()

test.__doc__ = "GEAR3 DF TRAN test"

if __name__ == '__main__':
    nt = NetlistTest('tran_gear3', er=1e-3, ea=1e-3)
    nt.setUp()
    nt.test()
