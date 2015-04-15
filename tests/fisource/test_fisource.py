from ahkab.testing import NetlistTest
from ahkab import options
options.plotting_show_plots = False

def test():
    nt = NetlistTest('fisource')
    nt.setUp()
    nt.test()
    nt.tearDown()

test.__doc__ = "CCCS test"

if __name__ == '__main__':
    nt = NetlistTest('fisource')
    nt.setUp()
    nt.test()
