from ahkab.testing import NetlistTest
from ahkab import options
options.plotting_show_plots = False

def test():
    nt = NetlistTest('bfpss1', er=1e-3, ea=1e-2)
    nt.setUp()
    nt.test()
    nt.tearDown()

test.__doc__ = "PSS BFPSS test 1"

if __name__ == '__main__':
    nt = NetlistTest('bfpss1')
    nt.setUp()
    nt.test()
