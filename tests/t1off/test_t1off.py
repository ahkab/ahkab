from ahkab.testing import NetlistTest
from ahkab import options
options.plotting_show_plots = False

def test():
    nt = NetlistTest('t1off', er=1e-3, ea=1e-5)
    nt.setUp()
    nt.test()
    nt.tearDown()

test.__doc__ = "SL: NMOS (OFF) with drain resistor."

if __name__ == '__main__':
    nt = NetlistTest('t1off')
    nt.setUp()
    nt.test()
