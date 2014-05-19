from ahkab.testing import NetlistTest
from ahkab import options
options.plotting_show_plots = False

def test():
    nt = NetlistTest('colpitts', er=1e-3, ea=1e-5)
    nt.setUp()
    nt.test()
    nt.tearDown()

test.__doc__ = "Transient simulation of a Colpitts oscillator"

if __name__ == '__main__':
    nt = NetlistTest('colpitts')
    nt.setUp()
    nt.test()
