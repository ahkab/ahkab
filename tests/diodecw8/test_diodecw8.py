from ahkab.testing import NetlistTest
from ahkab import options
options.plotting_show_plots = False

def myoptions():
    sim_opts = {}
    sim_opts.update({'vea':1e-9})
    sim_opts.update({'iea':1e-9})
    return sim_opts

def test():
    nt = NetlistTest('diodecw8', ea=1e-1, er=1e-3,
                     sim_opts=myoptions())
    nt.setUp()
    nt.test()
    if __name__ != '__main__':
        nt.tearDown()

test.__doc__ = "Diode CW multiplier x8 test"

if __name__ == '__main__':
    test()
