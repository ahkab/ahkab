from ahkab.testing import NetlistTest
from ahkab import options
options.plotting_show_plots = False

def test():
    nt = NetlistTest('downscaling_cm', ea=1e-6, er=1e-5)
    nt.setUp()
    nt.test()
    nt.tearDown()

test.__doc__ = "Downscaling current mirror (EKV) test"

if __name__ == '__main__':
    nt = NetlistTest('downscaling_cm', ea=1e-6, er=1e-5)
    nt.setUp()
    nt.test()
