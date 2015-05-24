import numpy as np
import ahkab
from ahkab import testing

cli = False

def test():
    """Test NL AC analysis (API)"""

    cir = ahkab.Circuit('CS amplifier')
    mys = ahkab.time_functions.sin(0, 1, 60)
    cir.add_vsource('vin', '1', '0', dc_value=3, ac_value=1)
    cir.add_vsource('vdd', '3', '0', dc_value=30)
    cir.add_resistor('Rd', '3', '2', 10e3)
    cir.add_capacitor('Cd', '3', '2', 40e-12)
    cir.add_resistor('Rs', '4', '0', 1e3)
    cir.add_capacitor('Cs', '4', '0', 4e-6)
    cir.add_model('ekv', 'ekv0', {'TYPE':'n', 'VTO':.4, 'KP':1e-2})
    cir.add_mos('m1', '2', '1', '4', '0', w=100e-6, l=1e-6, model_label='ekv0')
    print(cir)

    opa = ahkab.new_op(outfile='acnl', verbose=6)
    aca = ahkab.new_ac(1, 100e6, 10e3, outfile='acnl', verbose=6)
    r = ahkab.run(cir, [opa, aca])['ac']

    testbench = testing.APITest('acnl', cir, [opa, aca], skip_on_travis=False,
                                er=1e-3, ea=1e-5)
    testbench.setUp()
    testbench.test()

    if not cli:
        testbench.tearDown()

if __name__ == '__main__':
    cli = True
    test()
