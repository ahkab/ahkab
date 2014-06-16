from ahkab import options
options.plotting_show_plots = False
from ahkab.testing import NetlistTest

def test():
    nt = NetlistTest('pz1', er=1e-4, ea=1e-3)
    nt.setUp()
    nt.test()
    nt.tearDown()

test.__doc__ = "Pole Zero test 1"

if __name__ == '__main__':
    nt = NetlistTest('pz1', er=1e-4, ea=1e-3)
    nt.setUp()
    nt.test()
