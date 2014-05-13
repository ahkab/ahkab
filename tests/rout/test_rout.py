from ahkab.testing import NetlistTest
from ahkab import options
options.plotting_show_plots = False

def test():
    nt = NetlistTest('rout')
    nt.setUp()
    nt.test()
    nt.tearDown()

test.__doc__ = "Symbolic output resistance of a deg NMOS"

if __name__ == '__main__':
    nt = NetlistTest('rout')
    nt.setUp()
    nt.test()
