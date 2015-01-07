from ahkab.testing import NetlistTest
from ahkab import options
options.plotting_show_plots = False

def test():
    nt = NetlistTest('switch1')
    nt.setUp()
    nt.test()
    nt.tearDown()

test.__doc__ = "Switch test 1"

if __name__ == '__main__':
    nt = NetlistTest('switch1')
    nt.setUp()
    nt.test()
