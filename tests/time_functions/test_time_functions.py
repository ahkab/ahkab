import numpy as np
import ahkab
from ahkab import ahkab, circuit, printing, time_functions, testing

cli = False

def test():
    """Test pulse and sin API"""
    step = time_functions.pulse(v1=0, v2=1, td=500e-9, tr=1e-12, pw=1, tf=1e-12, per=2)
    damped_sin = time_functions.sin(vo=0, va=1, td=500e-9, freq=15e3, theta=5e3, phi=90.)
    exp = time_functions.exp(v1=.5, v2=-.05, td1=0, tau1=20e-6, td2=400e-6, tau2=20e-6)

    mycircuit = circuit.Circuit(title="Butterworth Example circuit", filename=None)

    gnd = mycircuit.get_ground_node()

    mycircuit.add_resistor(part_id="R1", n1="n1", n2="n2", value=600)
    mycircuit.add_inductor(part_id="L1", n1="n2", n2="n3", value=15.24e-3)
    mycircuit.add_capacitor(part_id="C1", n1="n3", n2=gnd, value=119.37e-9)
    mycircuit.add_inductor(part_id="L2", n1="n3", n2="n4", value=61.86e-3)
    mycircuit.add_capacitor(part_id="C2", n1="n4", n2=gnd, value=155.12e-9)
    mycircuit.add_resistor(part_id="R2", n1="n4", n2=gnd, value=1.2e3)

    mycircuit.add_vsource("V1", n1="n1", n2='n5', dc_value=3.3333, ac_value=.33333, function=step)
    mycircuit.add_vsource("V2", n1="n5", n2='n6', dc_value=3.3333, ac_value=.33333, function=damped_sin)
    mycircuit.add_vsource("V3", n1="n6", n2=gnd, dc_value=3.3333, ac_value=.33333, function=exp)

    op_analysis = ahkab.new_op(outfile='time_functions')
    ac_analysis = ahkab.new_ac(start=1e3, stop=1e5, points=100, outfile='time_functions')
    tran_analysis = ahkab.new_tran(tstart=0, tstop=1.2e-3, tstep=1e-6, x0=None, outfile='time_functions')

    testbench = testing.APITest('time_functions', mycircuit, 
                                [op_analysis, ac_analysis, tran_analysis],
                                skip_on_travis=True, er=1e-3, ea=1e-5)
    testbench.setUp()
    testbench.test()

    if cli:
        r = ahkab.run(mycircuit, an_list=[op_analysis, ac_analysis, tran_analysis])
        fig = plt.figure()
        plt.title(mycircuit.title + " - TRAN Simulation")
        plt.plot(r['tran']['T'], r['tran']['VN1'], label="Input voltage")
        plt.hold(True)
        plt.plot(r['tran']['T'], r['tran']['VN4'], label="output voltage")
        plt.legend()
        plt.hold(False)
        plt.grid(True)
        #plt.ylim([0,1.2])
        plt.ylabel('Step response')
        plt.xlabel('Time [s]')
        fig.savefig('tran_plot.png')

        fig = plt.figure()
        plt.subplot(211)
        plt.semilogx(r['ac']['w'], np.abs(r['ac']['Vn4']), 'o-')
        plt.ylabel('abs(V(n4)) [V]')
        plt.title(mycircuit.title + " - AC Simulation")
        plt.subplot(212)
        plt.grid(True)
        plt.semilogx(r['ac']['w'], np.angle(r['ac']['Vn4']), 'o-')
        plt.xlabel('Angular frequency [rad/s]')
        plt.ylabel('arg(V(n4)) [rad]')
        fig.savefig('ac_plot.png')
    else:
        testbench.tearDown()

if __name__ == '__main__':
    import pylab as plt
    cli = True
    test()
    plt.show()
