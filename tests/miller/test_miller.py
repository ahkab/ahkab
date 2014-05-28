from ahkab.testing import NetlistTest
from ahkab import options
options.plotting_show_plots = False

def test():
    nt = NetlistTest('miller')
    nt.setUp()
    nt.test()
    nt.tearDown()

test.__doc__ = "Symbolic Miller pole-splitting test"

if __name__ == '__main__':
    nt = NetlistTest('miller')
    nt.setUp()
    nt.test()
