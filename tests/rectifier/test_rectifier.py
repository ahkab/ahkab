import ahkab
from ahkab import circuit, printing, devices, testing

cli = False

def assemble():
    RS = 1
    cir = ahkab.circuit.Circuit("Rectifier", filename=None)
    sinf = ahkab.devices.sin(0, 220, 50.)
    cir.add_vsource('V1', 'inar', 'inbr', 1., 0., sinf)
    cir.add_resistor('R0a', 'inar', 'ina', 1.)
    cir.add_resistor('R0b', 'inbr', 'inb', .001)
    cir.add_inductor('L1a', 'ina', cir.gnd, 100)
    cir.add_inductor('L1b', cir.gnd, 'inb', 100)
    cir.add_inductor('L2a', 'int1', cir.gnd, .5)
    cir.add_inductor('L2b', cir.gnd, 'int2', .5)
    cir.add_inductor_coupling('K1a', 'L1a', 'L2a', .98)
    cir.add_inductor_coupling('K1b', 'L1b', 'L2b', .98)
    cir.add_model('diode', 'pnj', {'IS':1e-9, 'RS':RS})
    cir.add_diode('D1', 'int1', 'int4', 'pnj')
    cir.add_diode('D2', 'int3', 'int2', 'pnj')
    cir.add_diode('D3', 'int2', 'int4', 'pnj')
    cir.add_diode('D4', 'int3', 'int1', 'pnj')
#    cir.add_diode('D1', 'int1', 'int4f2', 'pnj')
#    cir.add_diode('D2', 'int3', 'int2f', 'pnj')
#    cir.add_diode('D3', 'int2', 'int4f', 'pnj')
#    cir.add_diode('D4', 'int3', 'int1f', 'pnj')
#    cir.add_resistor('R1D', 'int1f', 'int1', 1)
#    cir.add_resistor('R2D', 'int2f', 'int2', 1)
#    cir.add_resistor('R3D', 'int4f2', 'int4', 1)
#    cir.add_resistor('R4D', 'int4f', 'int4', 1)
    cir.add_capacitor('C1a', 'int3', 'int4', 100e-12)
    cir.add_capacitor('C1b', 'int3', 'int4', 10e-6)
    cir.add_resistor('R1', 'int3', 'int4', 1e4)
    cir.add_resistor('R1d', 'int1', 'int4', 1e5)
    cir.add_resistor('R2d', 'int3', 'int2', 1e5)
    cir.add_resistor('R3d', 'int2', 'int4', 1e5)
    cir.add_resistor('R4d', 'int3', 'int1', 1e5)

    if cli:
        printing.print_circuit(cir)

    return cir

def test():
    """Full wave rectifier test circuit"""

    cir = assemble()

    ## define analyses
    op1 = ahkab.new_op(outfile='rectifier')
    tran1 = ahkab.new_tran(0, 200e-3, 1e-4, outfile='rectifier', 
                           verbose=0+cli*6)

    # set the options
    sim_opts = {}
    sim_opts.update({'gmin':1e-7})
    sim_opts.update({'nl_voltages_lock':False})
    sim_opts.update({'nl_voltages_lock_factor':20})
    sim_opts.update({'iea':1e-1})
    sim_opts.update({'default_tran_method':'TRAP'})
    sim_opts.update({'hmin':1e-20})
    sim_opts.update({'transient_max_nr_iter':200})


    ## create a testbench
    testbench = testing.APITest('rectifier', cir, [op1, tran1],
                                skip_on_travis=True, sim_opts=sim_opts,
                                ea=1e-1, er=1.)

    ## setup and test
    testbench.setUp()
    testbench.test()

    ## this section is recommended. If something goes wrong, you may call the
    ## test from the cli and the plots to video in the following will allow
    ## for quick inspection
    if cli:
        ## re-run the test to grab the results
        cir = assemble()
        res = ahkab.run(cir, an_list=[op1, tran1])
        # print-out for good measure
        print "OP Results:"
        print list(res['op'].items())
        ## plot and save interesting data
        fig = plt.figure()
        plt.title(cir.title + " inputs")
        plt.plot(res['tran'].get_x(), res['tran']['VINA']-res['tran']['VINB'], label='Transf. input')
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
