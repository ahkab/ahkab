from ahkab.testing import NetlistTest

def test():
    nt = NetlistTest('rvtest2')
    nt.setUp()
    nt.test()
    nt.tearDown()

test.__doc__ = "1 resistor 1 voltage source simulation"

if __name__ == '__main__':
    nt = NetlistTest('rvtest2')
    nt.setUp()
    nt.test()
