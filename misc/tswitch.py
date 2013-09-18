import ahkab
import switch, circuit, printing
mycircuit = circuit.circuit(title="Test switch")
gnd = mycircuit.get_ground_node()
mycircuit.add_resistor("R1", "n1", gnd, 1e3)
mycircuit.add_resistor("R2", "n2", "n3", 50)
mycircuit.add_model("sw", "mysw", {'name':"mysw", 'VT':6, 'VH':1., 'RON':100., 'ROFF':None})
mycircuit.add_switch('SW1', 'n2', gnd, 'n1', gnd, ic=True, model_label='mysw')
mycircuit.add_isource('I1', 'n1', gnd, 1e-3)
mycircuit.add_vsource('V1', 'n3', gnd, 1.)
printing.print_circuit(mycircuit)
opa = {'type':'op', 'guess_label':None, 'guess':False}
ahkab.process_analysis([opa], mycircuit, '/tmp/out.data', 6)
