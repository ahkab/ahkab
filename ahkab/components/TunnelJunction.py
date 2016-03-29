from __future__ import (unicode_literals, absolute_import, division, print_function)
from .Component import Component
from math import exp, sqrt
import scipy.constants;
import scipy as np
from .. import utilities
from .. import options

class TunnelJunction(Component):
    """A diode element.

    **Parameters:**

    n1, n2 : string
        The diode anode and cathode.
    model : model instance
        The diode model providing the mathemathical modeling.
    ic : float
        The diode initial voltage condition for transient analysis
        (ie :math:`V_D = V(n_1) - V(n_2)` at :math:`t = 0`).

    The other are the physical parameters reported in the following table:

    +---------------+-------------------+-----------------------------------+
    | *Parameter*   | *Default value*   | *Description*                     |
    +===============+===================+===================================+
    | d             | 1 A               | Distance between contacts         |
    +---------------+-------------------+-----------------------------------+

    """

    def __init__(self, part_id, n1, n2, model, d=1, ic=None):
        #self.value = value
        self.n1 = n1
        self.n2 = n2
        self.part_id = part_id
        self.ic = ic
        self.coupling_devices = []
        self.is_nonlinear = True
        self.is_symbolic = True
        self.is_timedependent = False
        self._time_function = None
        self.ports = ((self.n1, self.n2),)
        
        self.d=d
        
        self.dc_guess = [0.425]
        class _dev_class(object):
            pass
        self.device = _dev_class()
        self.device.d = d
        self.device.last_vd = .425
        
        self.model = model

        if self.ic is not None:  # fixme
            print("(W): ic support in tunnel junctions is very experimental.")
            self.dc_guess = ic

    def i(self, op_index, ports_v, time=0):  # with gmin added
        '''It returns the current flowing into the element if the voltages specified in the voltages_vector are applied to its ports, at the time given.'''
        v = ports_v[0]
        i = self.model.get_i(self.model, v, self.device)
        return i

    def g(self, op_index, ports_v, port_index, time=0):
        '''returns the differential transconductance between the port at position port_index in the ports_vector (see point 2 above) and the element output current, when the operating point is specified by the voltages in the voltages_vector.'''
        if not port_index == 0:
            raise Exception("Attepted to evaluate a diode's gm on an unknown port.")
        gm = self.model.get_gm(self.model, op_index, utilities.tuplinator(ports_v), port_index, self.device)
        return gm
    
    def istamp(self, ports_v, time=0, reduced=True):
        """Get the current matrix

        A matrix corresponding to the current flowing in the element
        with the voltages applied as specified in the ``ports_v`` vector.

        **Parameters:**

        ports_v : list
            A list in the form: [voltage_across_port0, voltage_across_port1, ...]
        time: float
            the simulation time at which the evaluation is performed.
            It has no effect here. Set it to None during DC analysis.

        """
        v = ports_v[0]
        i = self.model.get_i(self.model, v, self.device)
        istamp = np.array((i, -i), dtype=np.float64)
        indices = ((self.n1 - 1*reduced, self.n2 - 1*reduced), (0, 0))
        if reduced:
            delete_i = [pos for pos, ix in enumerate(indices[0]) if ix == -1]
            istamp = np.delete(istamp, delete_i, axis=0)
            indices = tuple(zip(*[(ix, j) for ix, j in zip(*indices) if ix != -1]))
        return indices, istamp

    def gstamp(self, ports_v, time=0, reduced=True):
        """Returns the differential (trans)conductance wrt the port specified by port_index
        when the element has the voltages specified in ports_v across its ports,
        at (simulation) time.

        ports_v: a list in the form: [voltage_across_port0, voltage_across_port1, ...]
        port_index: an integer, 0 <= port_index < len(self.get_ports())
        time: the simulation time at which the evaluation is performed. Set it to
        None during DC analysis.
        """
        indices = ([self.n1 - 1]*2 + [self.n2 - 1]*2,
                   [self.n1 - 1, self.n2 - 1]*2)
        gm = self.model.get_gm(self.model, 0, utilities.tuplinator(ports_v), 0, self.device)
        if gm == 0:
            gm = options.gmin*2
        stamp = np.array(((gm, -gm),
                          (-gm, gm)), dtype=np.float64)
        if reduced:
            zap_rc = [pos for pos, i in enumerate(indices[1][:2]) if i == -1]
            stamp = np.delete(stamp, zap_rc, axis=0)
            stamp = np.delete(stamp, zap_rc, axis=1)
            indices = tuple(zip(*[(i, y) for i, y in zip(*indices) if (i != -1 and y != -1)]))
            stamp_flat = stamp.reshape(-1)
            stamp_folded = []
            indices_folded = []
            for ix, it in enumerate([(i, y) for i, y in zip(*indices)]):
                if it not in indices_folded:
                    indices_folded.append(it)
                    stamp_folded.append(stamp_flat[ix])
                else:
                    w = indices_folded.index(it)
                    stamp_folded[w] += stamp_flat[ix]
            indices = tuple(zip(*indices_folded))
            stamp = np.array(stamp_folded)
        return indices, stamp

    def __str__(self):
        rep = "%s d=%g" % (self.model.name, self.device.d)
        if self.ic is not None:
            rep = rep + " ic=" + str(self.ic)
        return rep

    def get_output_ports(self):
        return self.ports

    def get_drive_ports(self, op):
        if not op == 0:
            raise ValueError('Tunnel juinction %s has no output port %d' %
                             (self.part_id, op))
        return self.ports


    def get_op_info(self, ports_v_v):
        """Information regarding the Operating Point (OP)

        **Parameters:**

        ports_v : list of lists
            The parameter is to be set to ``[[v]]``, where ``v`` is the voltage
            applied to the diode terminals.

        **Returns:**

        op_keys : list of strings
            The labels corresponding to the numeric values in ``op_info``.
        op_info : list of floats
            The values corresponding to ``op_keys``.
        """
        vn1n2 = float(ports_v_v[0][0])
        idiode = self.i(0, (vn1n2,))
        gmdiode = self.g(0, (vn1n2,), 0)
        op_keys = ["Part ID", "V(n1-n2) [V]", "I(n1-n2) [A]", "P [W]",
                "gm [A/V]"]
        op_info = [self.part_id.upper(), vn1n2, idiode, vn1n2*idiode, gmdiode]
        return op_keys, op_info

        # for inductor
        # vn1n2 = float(ports_v[0][0])
        # energy = .5 * self.value * current**2
        # op_keys = ['Part ID', u"\u03d5(n1,n2) [Wb]", "I(n1->n2) [A]", "E [J]"]
        # op_info = [self.part_id.upper(), self.value*current, current, energy]
        # return op_keys, op_info
        
        

    def get_netlist_elem_line(self, nodes_dict):
        """A netlist line that, parsed, evaluates to the same instance

        **Parameters:**

        nodes_dict : dict
            The nodes dictionary of the circuit, so that the method
            can convert its internal node IDs to the corresponding
            external ones.

        **Returns:**

        ntlst_line : string
            The netlist line.
        """
        ext_n1, ext_n2 = nodes_dict[self.n1], nodes_dict[self.n2]
        ret = "%s %s %s %s" % (self.part_id, ext_n1, ext_n2, self.model.name)
        # append the optional part:
        # [<AREA=float> <T=float> <IC=float> <OFF=boolean>]
        ret += " d=%g" % self.device.d
        if self.ic is not None:
            ret += " IC=%g" % self.ic
        return ret
    

PHI_B_DEFAULT = 10  # eV



pi= scipy.constants.pi;
hbar=scipy.constants.hbar;
q = scipy.constants.elementary_charge;
eV= scipy.constants.eV;
m = scipy.constants.electron_mass;
#V = 1;
endiff=lambda V,PHI_B,fr: (PHI_B - fr*q*V) / 2
d_endiff=lambda V,fr: - fr*q / 2
expfr=lambda V,PHI_B,d,fr: exp(- (2 * sqrt(2 * m)) / hbar * sqrt(endiff(V,PHI_B,fr)) * d)
d_expfr=lambda V,PHI_B,d,fr: - (2 * sqrt(2 * m)) / hbar * d / sqrt(endiff(V,PHI_B,fr)) / 2 * d_endiff(V,fr)*exp(- (2 * sqrt(2 * m) / hbar * sqrt(endiff(V,PHI_B,fr)) * d))
term=lambda V,PHI_B,d,fr: endiff(V,PHI_B,fr)*expfr(V,PHI_B,d,fr)
d_term=lambda V,PHI_B,d,fr: d_endiff(V,fr)*expfr(V,PHI_B,d,fr) + d_expfr(V,PHI_B,d,fr)*endiff(V,PHI_B,fr)

class TunnelJunctionModel(object):
    """A tunnel junction model implementing the Simmons' formula

    Currently the capacitance modeling part is missing.

    The principal parameters are:

    +---------------+-------------------+-----------------------------------+
    | *Parameter*   | *Default value*   | *Description*                     |
    +===============+===================+===================================+
    | PHI_B         | 10 eV             | Potential barrier height          |
    +---------------+-------------------+-----------------------------------+
    """
    def __init__(self, name, PHI_B=None):
        self.name = name
        self.PHI_B = float(PHI_B) if PHI_B is not None else PHI_B_DEFAULT

    def print_model(self):
        strm = ".model tunnel_junction %s PHI_B=%g"
        print(strm % (self.name, self.PHI_B))

    @utilities.memoize
    def get_i(self, vext, dev):
        i = self._get_i(vext, dev)
        dev.last_vd = vext
        return i

    @utilities.memoize
    def get_gm(self, op_index, ports_v, port_index, dev):
        V = ports_v[port_index]
        return (q / (4 * pi ** 2 * hbar * (dev.d*1E-10) ** 2)) * (d_term(V,self.PHI_B*eV,(dev.d*1E-10),1) - d_term(V,self.PHI_B*eV,(dev.d*1E-10),- 1))

    def _safe_exp(self, x):
        return np.exp(x) if x < 70 else np.exp(70) + 10 * x

    def _get_i(self, V, dev):
        return (q / (4 * pi ** 2 * hbar * (dev.d*1E-10) ** 2))*(term(V,self.PHI_B*eV,(dev.d*1E-10),1) - term(V,self.PHI_B*eV,(dev.d*1E-10),- 1))

    def __str__(self):
        pass

