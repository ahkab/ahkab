import numpy as np
import ahkab
from ahkab import ahkab, circuit, printing, devices, testing

cli = False

def test():
    """Test VCCS (API)"""
    # The circuit is:
    # test for transconductors
    # va 1 2 type=idc idc=1m
    # r1 1 0 .5k
    # r2 2 0 .5k
    # g1 3 4 2 1 1e-3
    # r3 3 0 1k
    # r4 4 0 1k
    # .op
    # .symbolic

    mycircuit = circuit.Circuit(title="Test CCCS API", filename=None)
    gnd = mycircuit.get_ground_node()

    mycircuit.add_resistor(name="R1", n1="1", n2=gnd, value=500)
    mycircuit.add_resistor(name="R2", n1="2", n2=gnd, value=500)
    mycircuit.add_vccs('G1', n1='3', n2='4', sn1='2', sn2='1', value=1e-3)
    mycircuit.add_resistor(name="R3", n1="3", n2=gnd, value=1e3)
    mycircuit.add_resistor(name="R4", n1="4", n2=gnd, value=1e3)
    mycircuit.add_isource("IA", n1="1", n2='2', dc_value=1e-3)

    print(mycircuit)

    op_analysis = ahkab.new_op(outfile='gisource_api', verbose=6)
    symb_analysis = ahkab.new_symbolic(outfile='gisource_api', verbose=6)

    testbench = testing.APITest('gisource', mycircuit, 
                                [op_analysis, symb_analysis],
                                skip_on_travis=False, er=1e-3, ea=1e-5)
    testbench.setUp()
    testbench.test()

    if not cli:
        testbench.tearDown()

if __name__ == '__main__':
    cli = True
    test()
