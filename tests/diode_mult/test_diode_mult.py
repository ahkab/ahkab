from ahkab.testing import NetlistTest
from ahkab import options
options.plotting_show_plots = False

def test():
    nt = NetlistTest('diode_mult')
    nt.setUp()
    nt.test()
    nt.tearDown()

test.__doc__ = "Diode multiplier test"

if __name__ == '__main__':
    nt = NetlistTest('diode_mult')
    nt.setUp()
    nt.test()
