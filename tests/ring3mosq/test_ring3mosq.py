from ahkab.testing import NetlistTest
from ahkab import options
options.plotting_show_plots = False

def test():
    nt = NetlistTest('ring3mosq', er=1e-1, ea=1e-2)
    nt.setUp()
    nt.test()
    nt.tearDown()

test.__doc__ = "Transient simulation of a 3-stages ring oscillator (MOSQ)"

if __name__ == '__main__':
    nt = NetlistTest('ring3mosq')
    nt.setUp()
    nt.test()
