from ahkab.testing import NetlistTest
from ahkab import options
options.plotting_show_plots = False

def test():
    nt = NetlistTest('shooting1', ea=1e-3, er=1e-3)
    nt.setUp()
    nt.test()
    nt.tearDown()

test.__doc__ = "PSS SHOOTING test 1"

if __name__ == '__main__':
    nt = NetlistTest('shooting1', ea=1e-3, er=1e-3)
    nt.setUp()
    nt.test()
