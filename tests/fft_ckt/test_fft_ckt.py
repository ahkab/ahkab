from ahkab.testing import NetlistTest
from ahkab import options
options.plotting_show_plots = False

def test():
    nt = NetlistTest('fft_ckt')
    nt.setUp()
    nt.test()
    nt.tearDown()

test.__doc__ = ".FFT circuit test"

if __name__ == '__main__':
    nt = NetlistTest('fft_ckt')
    nt.setUp()
    nt.test()
