from __future__ import print_function
import numpy as np
import sympy
import ahkab
from ahkab import Circuit
from ahkab import printing
from ahkab import testing
from ahkab import symbolic


cli = False

def buildsvf(svf):
    gnd = svf.gnd
    ar = svf.add_resistor
    ac = svf.add_capacitor
    al = svf.add_inductor
    ao = lambda name, p, n: svf.add_vcvs(name, "U"+name[1:]+"o", gnd, p, n, 1e6)
    ar("Rin", gnd, "U1p", 4.7e3)
    ar("R00", "U1p", "U2o", 10e3)
    ar("R01", "in", "U1n", 10e3)
    ar("Rf1", "U1o", "U2n", 10e3)
    ar("R10", "U1o", "U1n", 10e3)
    ar("R11", "U1n", "U3o", 10e3)
    ar("Rf2", "U2o", "U3n", 10e3)
    ar("R02", "U2o", "out", 10e3)
    ac("C10", "U2o", "U2n", 15e-9)
    ac("C11", "U3o", "U3n", 15e-9)
    ao("E1", "U1p", "U1n")
    ao("E2", gnd, "U2n")
    ao("E3",  gnd, "U3n")

def test():
    """Test SVF Biquad"""
    mycircuit = Circuit(title="state variable filter")
    gnd = mycircuit.get_ground_node()
    buildsvf(mycircuit)
    mycircuit.add_vsource(part_id="V1", n1="in", n2=gnd, dc_value=5, ac_value=1)

    if cli:
        printing.print_circuit(mycircuit)

    subs = symbolic.parse_substitutions(('E2=E1', 'E3=E1', 'R01=R00', 'R02=R00',
                                         'R11=R00', 'R10=R00', 'C11=C10', 'Rf2=Rf1',
                                         'Rin=R00'))

    symbolic_sim = ahkab.new_symbolic(ac_enable=True, subs=subs, outfile='svf_biquad')
    ac_sim = ahkab.new_ac(start=0.1, stop=100e6, points=1000, x0=None, outfile='svf_biquad')

    testbench = testing.APITest('svf_biquad', mycircuit, [symbolic_sim, ac_sim],
                                skip_on_travis=True, er=1e-3, ea=1e-5)
    testbench.setUp()
    testbench.test()

    if cli:
        r = ahkab.run(mycircuit, [symbolic_sim, ac_sim])
        E = r['symbolic'][0].as_symbol('E1')
        out_hp = sympy.limit(r['symbolic'][0]['VU1o'], E, sympy.oo, '+')
        out_bp = sympy.limit(r['symbolic'][0]['VU2o'], E, sympy.oo, '+')
        out_lp = sympy.limit(r['symbolic'][0]['VU3o'], E, sympy.oo, '+')
        out_hp = out_hp.simplify()
        out_bp = out_bp.simplify()
        out_lp = out_lp.simplify()
        print("VU1o =", out_hp)
        print("VU2o =", out_bp)
        print("VU3o =", out_lp)

        w = sympy.Symbol('w')
        out_hp = out_hp.subs({r['symbolic'][0].as_symbol('RF1'):10e3,
                              r['symbolic'][0].as_symbol('C10'):15e-9,
                              r['symbolic'][0].as_symbol('V1'):1,
                              r['symbolic'][0].as_symbol('s'):1j*w,
                              })
        out_bp = out_bp.subs({r['symbolic'][0].as_symbol('RF1'):10e3,
                              r['symbolic'][0].as_symbol('C10'):15e-9,
                              r['symbolic'][0].as_symbol('V1'):1,
                              r['symbolic'][0].as_symbol('s'):1j*w,
                              })
        out_lp = out_lp.subs({r['symbolic'][0].as_symbol('RF1'):10e3,
                              r['symbolic'][0].as_symbol('C10'):15e-9,
                              r['symbolic'][0].as_symbol('V1'):1,
                              r['symbolic'][0].as_symbol('s'):1j*w,
                              })
        out_lp = sympy.lambdify((w,), out_lp, modules='numpy')
        out_bp = sympy.lambdify((w,), out_bp, modules='numpy')
        out_hp = sympy.lambdify((w,), out_hp, modules='numpy')
        ws = r['ac']['w'][::30]
        fig = plt.figure()
        plt.title(mycircuit.title)
        plt.subplot(211)
        plt.hold(True)
        plt.semilogx(r['ac']['w']/2./np.pi, 20*np.log10(np.abs(r['ac']['VU1o'])), label="HP output (AC)")
        plt.semilogx(r['ac']['w']/2./np.pi, 20*np.log10(np.abs(r['ac']['VU2o'])), label="BP output (AC)")
        plt.semilogx(r['ac']['w']/2./np.pi, 20*np.log10(np.abs(r['ac']['VU3o'])), label="LP output (AC)")
        plt.semilogx(ws/2./np.pi, 20*np.log10(np.abs(out_hp(ws))), 'v', label="HP output (SYMB)")
        plt.semilogx(ws/2./np.pi, 20*np.log10(np.abs(out_bp(ws))), 'v', label="BP output (SYMB)")
        plt.semilogx(ws/2./np.pi, 20*np.log10(np.abs(out_lp(ws))), 'v', label="LP output (SYMB)")
        plt.hold(False)
        plt.grid(True)
        plt.legend()
        plt.ylabel('Magnitude [dB]')
        plt.xlabel('Frequency [Hz]')
        plt.subplot(212)
        plt.hold(True)
        plt.semilogx(r['ac']['w']/2./np.pi, np.angle(r['ac']['VU1o']), label="HP output (AC)")
        plt.semilogx(r['ac']['w']/2./np.pi, np.angle(r['ac']['VU2o']), label="BP output (AC)")
        plt.semilogx(r['ac']['w']/2./np.pi, np.angle(r['ac']['VU3o']), label="LP output (AC)")
        plt.semilogx(ws/2./np.pi, np.angle(out_hp(ws)), 'v', label="HP output (SYMB)")
        plt.semilogx(ws/2./np.pi, np.angle(out_bp(ws)), 'v', label="BP output (SYMB)")
        plt.semilogx(ws/2./np.pi, np.angle(out_lp(ws)), 'v', label="LP output (SYMB)")
        plt.legend()
        plt.hold(False)
        plt.grid(True)
        #plt.ylim([0,1.2])
        plt.ylabel('Phase [rad]')
        plt.xlabel('Frequency [Hz]')
        fig.savefig('ac_plot.png')
    else:
        testbench.tearDown()

if __name__ == '__main__':
    import pylab as plt
    cli = True
    test()
    plt.show()
