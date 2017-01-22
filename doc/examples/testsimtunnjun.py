from ahkab import new_dc, run
from ahkab.circuit import Circuit
from ahkab.plotting import plot_results, show_plots # calls matplotlib for you
import scipy as np
from ahkab.components.TunnelJunction import *


# Define the circuit


tjm=TunnelJunctionModel("tjm")


cir = Circuit('Scanning Tunnel Microscope IV')
cir.add_model("diode","ddd",{"IS":0.000139*3*2.5,"N":160})

cir.add_vsource('V1','n1', cir.gnd, dc_value=1., ac_value=0.)

cir.add_diode('D1','n1','n2',"ddd")
cir.add_resistor('R1', 'n1', 'n3', 1.)
tj=TunnelJunction('tj1',cir.add_node('n2'),cir.add_node(cir.gnd),tjm,d=10)
cir.append(tj)

dc1 = new_dc(-2.1,3, 1e2, source="V1")
res = run(cir, dc1)
plot_results('Scanning Tunnel Microscope IV', [('I(V1)','')], res['dc'])

show_plots()
