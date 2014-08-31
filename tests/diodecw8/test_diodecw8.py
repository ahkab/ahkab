from ahkab.testing import NetlistTest
from ahkab import options
options.plotting_show_plots = False

def test():
    nt = NetlistTest('diodecw8')
    nt.setUp()
    nt.test()
    nt.tearDown()

test.__doc__ = "Diode CW multiplier x8 test"

if __name__ == '__main__':
    nt = NetlistTest('diodecw8')
    nt.setUp()
    nt.test()
