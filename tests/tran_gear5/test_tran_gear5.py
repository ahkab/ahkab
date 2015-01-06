from __future__ import unicode_literals, print_function
from ahkab.testing import NetlistTest
from ahkab import options
options.plotting_show_plots = False

def test():
    nt = NetlistTest('tran_gear5', er=1e-2, ea=1e-2)
    nt.setUp()
    nt.test()
    nt.tearDown()

test.__doc__ = "GEAR5 DF TRAN test"

if __name__ == '__main__':
    nt = NetlistTest('tran_gear5', er=1e-2, ea=1e-3)
    nt.setUp()
    nt.test()
