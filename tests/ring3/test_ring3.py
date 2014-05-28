from ahkab.testing import NetlistTest
from ahkab import options
options.plotting_show_plots = False

def test():
    nt = NetlistTest('ring3', er=1e-1, ea=1e-2)
    nt.setUp()
    nt.test()
    nt.tearDown()

test.__doc__ = "Transient simulation of a 3-stages ring oscillator"

if __name__ == '__main__':
    nt = NetlistTest('ring3')
    nt.setUp()
    nt.test()
