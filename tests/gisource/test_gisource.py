from ahkab.testing import NetlistTest
from ahkab import options
options.plotting_show_plots = False

def test():
    nt = NetlistTest('gisource')
    nt.setUp()
    nt.test()
    nt.tearDown()

test.__doc__ = "VCCS test (netlist)"

if __name__ == '__main__':
    nt = NetlistTest('gisource')
    nt.setUp()
    nt.test()
