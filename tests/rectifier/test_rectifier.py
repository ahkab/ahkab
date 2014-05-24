import ahkab
from ahkab import circuit, printing, devices, testing

cli = False

def test():
    """Full wave rectifier test circuit"""

    cir = ahkab.circuit.Circuit("Rectifier", filename=None)
    sinf = ahkab.devices.sin(0, 220, 50.)
    cir.add_vsource('V1', 'inr', cir.gnd, 1., 0., sinf)
    cir.add_resistor('R0', 'inr', 'in', .001)
    cir.add_inductor('L1', 'in', cir.gnd, 100e-6)
    cir.add_inductor('L2', 'int1', 'int2', 100e-6)
    cir.add_inductor_coupling('K1', 'L1', 'L2', .04)
    cir.add_model('diode', 'pnj', {'IS':1e-10})
    cir.add_diode('D1', 'int1', 'int4', 'pnj')
    cir.add_diode('D2', 'int3', 'int2', 'pnj')
    cir.add_diode('D3', 'int2', 'int4', 'pnj')
    cir.add_diode('D4', 'int3', 'int1', 'pnj')
    cir.add_capacitor('C1a', 'int3', 'int4', 100e-12)
    cir.add_capacitor('C1b', 'int3', 'int4', 100e-6)
    cir.add_resistor('R1', 'int3', 'int4', 1000)


    if cli:
        printing.print_circuit(cir)

    ## define analyses
    op1 = ahkab.new_op(outfile='rectifier')
    tran1 = ahkab.new_tran(0, 200e-3, 1e-3, outfile='rectifier')

    ## create a testbench
    testbench = testing.APITest('rectifier', cir, [op1, tran1],
                                skip_on_travis=False, ea=1e-2, er=1e-2)

    ahkab.options.default_tran_method = "TRAP"
    ahkab.options.hmin = 1e-20
    ahkab.options.transient_max_nr_iter = 100

    ## setup and test
    testbench.setUp()
    testbench.test()

    ## this section is recommended. If something goes wrong, you may call the
    ## test from the cli and the plots to video in the following will allow
    ## for quick inspection
    if cli:
        ## re-run the test to grab the results
        res = ahkab.run(cir, an_list=[op1, tran1])
        # print-out for good measure
        print "OP Results:"
        print list(res['op'].items())
        ## plot and save interesting data
        fig = plt.figure()
        plt.title(cir.title + " inputs")
        plt.plot(res['tran'].get_x(), res['tran']['VIN'], label='Input voltage')
        plt.hold(True)
        plt.plot(res['tran'].get_x(), res['tran']['vint1'], label='Transformer output #1')
        plt.plot(res['tran'].get_x(), res['tran']['vint2'], label='Transformer output #2')
        plt.hold(False)
        plt.grid(True)
        plt.legend()
        plt.ylabel('Voltage [V]')
        plt.xlabel('Time [s]')
        fig.savefig('rectf1_plot.png')
        fig = plt.figure()
        plt.title(cir.title + " outputs")
        plt.plot(res['tran'].get_x(), res['tran']['vint4']-res['tran']['vint3'],
                 label="output voltage")
        plt.legend()
        plt.grid(True)
        plt.ylabel('Voltage [V]')
        plt.xlabel('Time [s]')
        fig.savefig('rectf2_plot.png')

    else:
        testbench.tearDown()

if __name__ == '__main__':
    import pylab as plt
    cli = True
    test()
    plt.show()
