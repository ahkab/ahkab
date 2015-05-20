from ahkab.testing import NetlistTest
from ahkab import options
options.plotting_show_plots = False

def test():
    nt = NetlistTest('fourier_ckt')
    nt.setUp()
    nt.test()
    nt.tearDown()

test.__doc__ = ".FOUR circuit test"

if __name__ == '__main__':
    nt = NetlistTest('fourier_ckt')
    nt.setUp()
    nt.test()
