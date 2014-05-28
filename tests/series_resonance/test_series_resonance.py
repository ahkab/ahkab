from ahkab.testing import NetlistTest
from ahkab import options
options.plotting_show_plots = False

def test():
    nt = NetlistTest('series_resonance')
    nt.setUp()
    nt.test()
    nt.tearDown()

test.__doc__ = "AC series resonance test"

if __name__ == '__main__':
    nt = NetlistTest('series_resonance')
    nt.setUp()
    nt.test()
