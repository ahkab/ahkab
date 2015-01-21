import numpy as np
import ahkab
from ahkab import ahkab, circuit, printing, devices, testing

cli = False

def test():
    """Test CCCS API"""
    # The circuit is:
    # test for transconductors
    # va 1 2 type=vdc vdc=.1
    # r1 1 0 .5k
    # r2 2 0 .5k
    # f1 3 4 va 5
    # r3 3 0 1k
    # r4 4 0 1k
    # .op
    # .symbolic

    mycircuit = circuit.Circuit(title="Test CCCS API", filename=None)
    gnd = mycircuit.get_ground_node()

    mycircuit.add_resistor(part_id="R1", n1="1", n2=gnd, value=500)
    mycircuit.add_resistor(part_id="R2", n1="2", n2=gnd, value=500)
    mycircuit.add_cccs('F1', n1='3', n2='4', source_id='VA', value=5)
    mycircuit.add_resistor(part_id="R3", n1="3", n2=gnd, value=1e3)
    mycircuit.add_resistor(part_id="R4", n1="4", n2=gnd, value=1e3)
    mycircuit.add_vsource("VA", n1="1", n2='2', dc_value=0.1)

    print(mycircuit)

    op_analysis = ahkab.new_op(outfile='fisource_api', verbose=6)
    symb_analysis = ahkab.new_symbolic(outfile='fisource_api', verbose=6)

    testbench = testing.APITest('fisource', mycircuit, 
                                [op_analysis, symb_analysis],
                                skip_on_travis=False, er=1e-3, ea=1e-5)
    testbench.setUp()
    testbench.test()

    if not cli:
        testbench.tearDown()

if __name__ == '__main__':
    cli = True
    test()
