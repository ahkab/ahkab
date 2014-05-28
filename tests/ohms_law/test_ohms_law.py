from ahkab.testing import NetlistTest
from ahkab import options
options.plotting_show_plots = False

def test():
    nt = NetlistTest('ohms_law')
    nt.setUp()
    nt.test()
    nt.tearDown()

test.__doc__ = "Ohm's law OP test"

if __name__ == '__main__':
    nt = NetlistTest('ohms_law')
    nt.setUp()
    nt.test()
