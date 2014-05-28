from ahkab.testing import NetlistTest
from ahkab import options
options.plotting_show_plots = False

def test():
    nt = NetlistTest('v2')
    nt.setUp()
    nt.test()
    nt.tearDown()

test.__doc__ = "OP double vsource test"

if __name__ == '__main__':
    nt = NetlistTest('v2')
    nt.setUp()
    nt.test()
