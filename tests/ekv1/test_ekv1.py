from ahkab.testing import NetlistTest
from ahkab import options
options.plotting_show_plots = False

def test():
    nt = NetlistTest('ekv1')
    nt.setUp()
    nt.test()
    nt.tearDown()

test.__doc__ = "EKV NMOS DC sweep"

if __name__ == '__main__':
    nt = NetlistTest('ekv1')
    nt.setUp()
    nt.test()
