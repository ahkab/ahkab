from ahkab.testing import NetlistTest
from ahkab import options

def test():
    nt = NetlistTest('symbro')
    nt.setUp()
    nt.test()
    nt.tearDown()

test.__doc__ = "OP double vsource test"

if __name__ == '__main__':
    nt = NetlistTest('symbro')
    nt.setUp()
    nt.test()
