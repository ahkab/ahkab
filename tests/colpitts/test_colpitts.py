from ahkab.testing import NetlistTest
from ahkab import options
options.plotting_show_plots = False

def myoptions():
        # optionally, set non-standard options
        sim_opts = {}
        sim_opts.update({'vea':1e-6})
        sim_opts.update({'iea':1e-6})
        #sim_opts.update({'transient_max_nr_iter':200})
        return sim_opts

def test():
    nt = NetlistTest('colpitts', er=1e-3, ea=5e-4, sim_opts=myoptions())
    nt.setUp()
    nt.test()
    nt.tearDown()

test.__doc__ = "Transient simulation of a Colpitts oscillator"

if __name__ == '__main__':
    nt = NetlistTest('colpitts', sim_opts=myoptions())
    nt.setUp()
    nt.test()
