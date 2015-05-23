import numpy as np
import ahkab
from ahkab import ahkab, circuit, printing, devices, testing

cli = False

def test():
    """Test CCCS API"""
    # The circuit is:
    #test for transresitances
    #va 1 2 type=vdc vdc=.1 vac=1
    #r1 1 0 .5k
    #r2 2 0 .5k
    #h1 3 4 va 5000
    #r3 3 0 1k
    #r4 4 5 1k
    #l1 5 0 10u
    #c1 5 0 10u
    #.op
    #.ac start=50k stop=5e5 nsteps=1000
    #.symbolic
    #.plot ac |v(5)|

    mycircuit = circuit.Circuit(title="Test CCVS API", filename=None)
    gnd = mycircuit.get_ground_node()

    mycircuit.add_resistor(part_id="R1", n1="1", n2=gnd, value=500)
    mycircuit.add_resistor(part_id="R2", n1="2", n2=gnd, value=500)
    mycircuit.add_vsource("VA", n1="1", n2='2', dc_value=0.1, ac_value=1.)
    mycircuit.add_ccvs('H1', n1='3', n2='4', source_id='VA', value=5000)
    mycircuit.add_resistor(part_id="R3", n1="3", n2=gnd, value=1e3)
    mycircuit.add_resistor(part_id="R4", n1="4", n2="5", value=1e3)
    mycircuit.add_inductor(part_id="L1", n1="5", n2=gnd, value=10e-6)
    mycircuit.add_capacitor(part_id="C1", n1="5", n2=gnd, value=10e-6)

    print(mycircuit)

    op_analysis = ahkab.new_op(outfile='hvsource_api', verbose=6)
    symb_analysis = ahkab.new_symbolic(outfile='hvsource_api', verbose=6)
    ac_analysis = ahkab.new_ac(outfile='hvsource_api', start=7957.747,
                               stop=79577.471, points=1000, verbose=6)

    testbench = testing.APITest('hvsource', mycircuit, 
                                [op_analysis, symb_analysis, ac_analysis],
                                skip_on_travis=False, er=1e-3, ea=1e-5)
    testbench.setUp()
    testbench.test()

    if not cli:
        testbench.tearDown()

if __name__ == '__main__':
    cli = True
    test()
