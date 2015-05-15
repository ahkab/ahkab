from ahkab.testing import NetlistTest
from ahkab import options
options.plotting_show_plots = False

def test():
    nt = NetlistTest('sffmckt')
    nt.setUp()
    nt.test()
    nt.tearDown()

test.__doc__ = "SFFM circuit test"

if __name__ == '__main__':
    nt = NetlistTest('sffmckt')
    nt.setUp()
    nt.test()
