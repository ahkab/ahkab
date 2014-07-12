from ahkab.testing import NetlistTest
from ahkab import options
options.plotting_show_plots = False

def test():
    nt = NetlistTest('t3n', er=1e-3, ea=1e-5)
    nt.setUp()
    nt.test()
    nt.tearDown()

test.__doc__ = "SL: NMOS (ON) PMOS (OFF) CMOS inverter."

if __name__ == '__main__':
    nt = NetlistTest('t3n')
    nt.setUp()
    nt.test()
