from __future__ import unicode_literals, print_function
from ahkab.testing import NetlistTest
from ahkab import options
options.plotting_show_plots = False

def test():
    nt = NetlistTest('tran_ie')
    nt.setUp()
    nt.test()
    nt.tearDown()

test.__doc__ = "IMPLICIT_EULER DF TRAN test"

if __name__ == '__main__':
    nt = NetlistTest('tran_ie')
    nt.setUp()
    nt.test()
