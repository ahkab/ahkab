import numpy as np
import ahkab
from ahkab import ahkab, circuit, time_functions

mycircuit = circuit.Circuit(title="Butterworth Example circuit", filename=None)

gnd = mycircuit.get_ground_node()

mycircuit.add_resistor("R1", n1="n1", n2="n2", value=600)
mycircuit.add_inductor("L1", n1="n2", n2="n3", value=15.24e-3)
mycircuit.add_capacitor("C1", n1="n3", n2=gnd, value=119.37e-9)
mycircuit.add_inductor("L2", n1="n3", n2="n4", value=61.86e-3)
mycircuit.add_capacitor("C2", n1="n4", n2=gnd, value=155.12e-9)
mycircuit.add_resistor("R2", n1="n4", n2=gnd, value=1.2e3)

voltage_step = time_functions.pulse(v1=0, v2=1, td=500e-9, tr=1e-12, pw=1, tf=1e-12, per=2)
mycircuit.add_vsource("V1", n1="n1", n2=gnd, dc_value=5, ac_value=1, function=voltage_step)

print mycircuit

op_analysis = ahkab.new_op()
ac_analysis = ahkab.new_ac(start=1e3, stop=1e5, points=100)
tran_analysis = ahkab.new_tran(tstart=0, tstop=1.2e-3, tstep=1e-6, x0=None)

r = ahkab.run(mycircuit, an_list=[op_analysis, ac_analysis, tran_analysis])

import pylab

fig = pylab.figure()
pylab.title(mycircuit.title + " - TRAN Simulation")
pylab.plot(r['tran']['T'], r['tran']['VN1'], label="Input voltage")
pylab.hold(True)
pylab.plot(r['tran']['T'], r['tran']['VN4'], label="output voltage")
pylab.legend()
pylab.hold(False)
pylab.grid(True)
pylab.ylim([0,1.2])
pylab.ylabel('Step response')
pylab.xlabel('Time [s]')
fig.savefig('tran_plot.png')

fig = pylab.figure()
pylab.subplot(211)
pylab.semilogx(r['ac']['f'], np.abs(r['ac']['Vn4']), 'o-')
pylab.ylabel('abs(V(n4)) [V]')
pylab.title(mycircuit.title + " - AC Simulation")
pylab.subplot(212)
pylab.grid(True)
pylab.semilogx(r['ac']['f'], np.angle(r['ac']['Vn4']), 'o-')
pylab.xlabel('Frequency [Hz]')
pylab.ylabel('arg(V(n4)) [rad]')
fig.savefig('ac_plot.png')
pylab.show()
