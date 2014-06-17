# -*- coding: iso-8859-1 -*-
# symbolic.py
# Symbolic simulation module
# Copyright 2010-2013 Giuseppe Venturini

# This file is part of the ahkab simulator.
#
# Ahkab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 2 of the License.
#
# Ahkab is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License v2
# along with ahkab.  If not, see <http://www.gnu.org/licenses/>.

"""
This module offers the functions needed to perform a symbolic simulation,
AC or DC.

The principal method is solve(), which carries out the symbolic circuit solution.
"""

import sympy
from sympy.matrices import zeros as smzeros

from . import circuit
from . import devices
from . import ekv
from . import mosq
from . import diode
from . import printing
from . import results
from . import options

specs = {'symbolic': {'tokens': ({
                                  'label': 'tf',
                                  'pos': None,
                                  'type': str,
                                  'needed': False,
                                  'dest': 'source',
                                  'default': None
                                 },
                                 {
                                  'label': 'ac',
                                  'pos': None,
                                  'type': bool,
                                  'needed': False,
                                  'dest': 'ac_enable',
                                  'default': True
                                 },
                                 {
                                  'label': 'r0s',
                                  'pos': None,
                                  'type': bool,
                                  'needed': False,
                                  'dest': 'r0s',
                                  'default': False
                                 }
                                )
                     }
        }

# the s variable
s = sympy.Symbol('s', complex=True)

enabled_assumptions = {'real':False, 'positive':False, 'complex':True}

def symbolic_analysis(circ, source=None, ac_enable=True, r0s=False, subs=None, outfile=None, verbose=3):
    """Attempt a symbolic solution of the circuit.

    **Parameters:**

    circ : circuit instance
        the circuit instance to be simulated.

    source : string, optional
        the ``part_id`` of the source to be used as input for the transfer
           function. If ``None``, no transfer function is evaluated.

    ac_enable : bool, optional
        take frequency dependency into consideration (default: True).

    r0s : bool, optional
        take transistors' output impedance into consideration (default: False)

    subs: dict, optional
        a dictionary of ``sympy.Symbol`` to be substituted. It makes solving the circuit
        easier. Eg. ``subs={R1:R2}`` - replace R1 with R2. It can be generated with
        :func:`parse_substitutions()`

    outfile : string, optional
        output filename - ``'stdout'`` means print to stdout, the default.

    verbose: int, optional
        verbosity level 0 (silent) to 6 (painful).

    **Returns:** 

    sol : dict
        a dictionary with the solutions.

    """
    if subs is None:
        subs = {}  # no subs by default

    if not ac_enable:
        printing.print_info_line(
            ("Starting symbolic DC analysis...", 1), verbose)
    else:
        printing.print_info_line(
            ("Starting symbolic AC analysis...", 1), verbose)

    printing.print_info_line(
        ("Building symbolic MNA, N and x...", 3), verbose, print_nl=False)
    mna, N, subs_g = generate_mna_and_N(
        circ, opts={'r0s': r0s}, ac=ac_enable, verbose=verbose)
    x = get_variables(circ)
    mna = mna[1:, 1:]
    N = N[1:, :]
    printing.print_info_line((" done.", 3), verbose)

    printing.print_info_line(
        ("Performing variable substitutions...", 5), verbose)
    mna, N = apply_substitutions(mna, N, subs)

    printing.print_info_line(("MNA matrix (reduced):", 5), verbose)
    printing.print_info_line((sympy.sstr(mna), 5), verbose)
    printing.print_info_line(("N matrix (reduced):", 5), verbose)
    printing.print_info_line((sympy.sstr(N), 5), verbose)

    printing.print_info_line(("Building equations...", 3), verbose)
    eq = []
    for i in _to_real_list(mna * x + N):
        eq.append(sympy.Eq(i, 0))

    x = _to_real_list(x)

    if verbose > 3:
        printing.print_symbolic_equations(eq)
    printing.print_info_line(("To be solved for:", 4), verbose)
    printing.print_info_line((str(x), 4), verbose)
    printing.print_info_line(("Solving...", 1), verbose)

    sol = sympy.solve(
            eq, x, manual=options.symb_sympy_manual_solver, simplify=True)

    for ks in sol.keys():
        sol.update({ks: sol[ks].subs(subs_g)})

    # sol = sol_to_dict(sol, x)

    if sol == {}:
        printing.print_warning("No solutions. Check the netlist.")
    else:
        printing.print_info_line(("Success!", 2), verbose)
        printing.print_info_line(("Results:", 1), verbose)
        if options.cli:
            printing.print_symbolic_results(sol)

    if source is not None:
        src = _symbol_factory(source.upper())
        printing.print_info_line(("Calculating small-signal symbolic transfer functions (%s))..." %
                                 (str(src),), 2), verbose, print_nl=False)
        tfs = calculate_gains(sol, src)
        printing.print_info_line(("done.", 2), verbose)
        printing.print_info_line(
            ("Small-signal symbolic transfer functions:", 1), verbose)
        if options.cli:
            printing.print_symbolic_transfer_functions(tfs)
    else:
        tfs = None

    # convert to a results instance
    sol = results.symbolic_solution(sol, subs, circ, outfile)
    if tfs:
        if outfile and outfile != 'stdout':
            outfile += ".tfs"
        tfs = results.symbolic_solution(tfs, subs, circ, outfile, tf=True)
    return sol, tfs


def calculate_gains(sol, xin, optimize=True):
    gains = {}
    for key, value in sol.iteritems():
        tf = {}
        gain = sympy.together(value.diff(xin)) if optimize else value.diff(xin)
        (ps, zs) = get_roots(gain)
        tf.update({'gain': gain})
        tf.update({'gain0': gain.subs(s, 0)})
        tf.update({'poles': ps})
        tf.update({'zeros': zs})
        gains.update({"%s/%s" % (str(key), str(xin)): tf})
    return gains


def sol_to_dict(sol, x, optimize=True):
    ret = {}
    for index in range(x.shape[0]):
        sol_current = sympy.together(sol[index]) if optimize else sol[index]
        ret.update({str(x[index]): sol_current})
    return ret


def apply_substitutions(mna, N, subs):
    mna = mna.subs(subs)
    N = N.subs(subs)
    return (mna, N)


def get_variables(circ):
    """Returns a sympy matrix with the circuit variables to be solved for.
    """
    nv_1 = len(circ.nodes_dict) - \
        1  # numero di soluzioni di tensione (al netto del ref)

    # descrizioni dei componenti non definibili in tensione
    idescr = [elem.part_id.upper()
              for elem in circ if circuit.is_elem_voltage_defined(elem)]

    mna_size = nv_1 + len(idescr)
    x = smzeros((mna_size, 1))

    for i in range(mna_size):
        if i < nv_1:
            x[i, 0] = _symbol_factory("V" + str(circ.nodes_dict[i + 1]))
        else:
            x[i, 0] = _symbol_factory("I[" + idescr[i - nv_1] + "]")
    return x


def _to_real_list(M):
    """
    M.tolist() returns a list of lists, even when the symb matrix is really just a vector.
    we want a list of symbols! This fixes that.

    mylist[k] = mymat.tolist[k][0]

    M: a sympy matrix with only one column

    Returns: a list.
    """
    fakelist = M.tolist()
    reallist = []
    for elem in fakelist:
        reallist.append(elem[0])
    return reallist


def generate_mna_and_N(circ, opts, ac=False, verbose=3):
    """Generates a symbolic Modified Nodal Analysis matrix and N vector.
    """
    #   print options
    n_of_nodes = len(circ.nodes_dict)
    mna = smzeros(n_of_nodes)
    N = smzeros((n_of_nodes, 1))
    subs_g = {}

    for elem in circ:
        if isinstance(elem, devices.Resistor):
            # we use conductances instead of 1/R because there is a significant
            # overhead handling many 1/R terms in sympy.
            if elem.is_symbolic:
                R = _symbol_factory(
                    elem.part_id.upper(), real=True, positive=True)
                G = _symbol_factory('G' + elem.part_id[1:], real=True, positive=True)
                # but we keep track of which is which and substitute back after
                # solving.
                subs_g.update({G: 1 / R})
            else:
                R = elem.value
                G = 1.0 / R
            mna[elem.n1, elem.n1] = mna[elem.n1, elem.n1] + G
            mna[elem.n1, elem.n2] = mna[elem.n1, elem.n2] - G
            mna[elem.n2, elem.n1] = mna[elem.n2, elem.n1] - G
            mna[elem.n2, elem.n2] = mna[elem.n2, elem.n2] + G
        elif isinstance(elem, devices.Capacitor):
            if ac:
                if elem.is_symbolic:
                    capa = _symbol_factory(
                        elem.part_id.upper(), real=True, positive=True)
                else:
                    capa = elem.value
                mna[elem.n1, elem.n1] = mna[elem.n1, elem.n1] + s * capa
                mna[elem.n1, elem.n2] = mna[elem.n1, elem.n2] - s * capa
                mna[elem.n2, elem.n2] = mna[elem.n2, elem.n2] + s * capa
                mna[elem.n2, elem.n1] = mna[elem.n2, elem.n1] - s * capa
            else:
                pass
        elif isinstance(elem, devices.Inductor):
            pass
        elif isinstance(elem, devices.GISource):
            if elem.is_symbolic:
                alpha = _symbol_factory(elem.part_id.upper(), real=True)
            else:
                alpha = elem.value
            mna[elem.n1, elem.sn1] = mna[elem.n1, elem.sn1] + alpha
            mna[elem.n1, elem.sn2] = mna[elem.n1, elem.sn2] - alpha
            mna[elem.n2, elem.sn1] = mna[elem.n2, elem.sn1] - alpha
            mna[elem.n2, elem.sn2] = mna[elem.n2, elem.sn2] + alpha
        elif isinstance(elem, devices.ISource):
            if elem.is_symbolic:
                IDC = _symbol_factory(elem.part_id.upper())
            else:
                IDC = elem.dc_value
            N[elem.n1, 0] = N[elem.n1, 0] + IDC
            N[elem.n2, 0] = N[elem.n2, 0] - IDC
        elif isinstance(elem, mosq.mosq_device) or isinstance(elem, ekv.ekv_device):
            gm = _symbol_factory('gm_' + elem.part_id, real=True, positive=True)
            mna[elem.n1, elem.ng] = mna[elem.n1, elem.ng] + gm
            mna[elem.n1, elem.n2] = mna[elem.n1, elem.n2] - gm
            mna[elem.n2, elem.ng] = mna[elem.n2, elem.ng] - gm
            mna[elem.n2, elem.n2] = mna[elem.n2, elem.n2] + gm
            if opts['r0s']:
                r0 = _symbol_factory(
                    'r0_' + elem.part_id, real=True, positive=True)
                mna[elem.n1, elem.n1] = mna[elem.n1, elem.n1] + 1 / r0
                mna[elem.n1, elem.n2] = mna[elem.n1, elem.n2] - 1 / r0
                mna[elem.n2, elem.n1] = mna[elem.n2, elem.n1] - 1 / r0
                mna[elem.n2, elem.n2] = mna[elem.n2, elem.n2] + 1 / r0
        elif isinstance(elem, diode.diode):
            gd = _symbol_factory("g" + elem.part_id, positive=True)
            mna[elem.n1, elem.n1] = mna[elem.n1, elem.n1] + gd
            mna[elem.n1, elem.n2] = mna[elem.n1, elem.n2] - gd
            mna[elem.n2, elem.n1] = mna[elem.n2, elem.n1] - gd
            mna[elem.n2, elem.n2] = mna[elem.n2, elem.n2] + gd
        elif isinstance(elem, devices.InductorCoupling):
            pass
            # this is taken care of within the inductors
        elif circuit.is_elem_voltage_defined(elem):
            pass
            # we'll add its lines afterwards
        elif verbose:
            printing.print_warning(
                "Skipped elem %s: not implemented." % (elem.part_id.upper(),))

    pre_vde = mna.shape[0]
    for elem in circ:
        if circuit.is_elem_voltage_defined(elem):
            index = mna.shape[0]  # get_matrix_size(mna)[0]
            mna = expand_matrix(mna, add_a_row=True, add_a_col=True)
            N = expand_matrix(N, add_a_row=True, add_a_col=False)
            # KCL
            mna[elem.n1, index] = +1
            mna[elem.n2, index] = -1
            # KVL
            mna[index, elem.n1] = +1
            mna[index, elem.n2] = -1
            if isinstance(elem, devices.VSource):
                if elem.is_symbolic:
                    VDC = _symbol_factory(elem.part_id.upper())
                else:
                    VDC = elem.dc_value
                N[index, 0] = -VDC
            elif isinstance(elem, devices.EVSource):
                if elem.is_symbolic:
                    alpha = _symbol_factory(elem.part_id.upper(), real=True)
                else:
                    alpha = elem.alpha
                mna[index, elem.sn1] = -alpha
                mna[index, elem.sn2] = +alpha
            elif isinstance(elem, devices.Inductor):
                if ac:
                    if elem.is_symbolic:
                        L = _symbol_factory(
                            elem.part_id.upper(), real=True, positive=True)
                    else:
                        L = elem.L
                    mna[index, index] = -s * L
                else:
                    pass
                    # already so: commented out
                    # N[index,0] = 0
            elif isinstance(elem, devices.HVSource):
                printing.print_warning(
                    "symbolic.py: BUG - hvsources are not implemented yet.")
                sys.exit(33)

    for elem in circ:
        if circuit.is_elem_voltage_defined(elem):
            if isinstance(elem, devices.Inductor):
                if ac:
                    # find its index to know which column corresponds to its
                    # current
                    this_index = circ.find_vde_index(elem.part_id, verbose=0)
                    for cd in elem.coupling_devices:
                        if cd.is_symbolic:
                            M = _symbol_factory(
                                cd.part_id, real=True, positive=True)
                        else:
                            M = cd.K
                        # get `part_id` of the other inductor (eg. "L32")
                        other_id_wdescr = cd.get_other_inductor(elem.part_id)
                        # find its index to know which column corresponds to
                        # its current
                        other_index = circ.find_vde_index(
                            other_id_wdescr, verbose=0)
                        # add the term.
                        mna[pre_vde + this_index,
                            pre_vde + other_index] += -s * M
                else:
                    pass

    # all done
    return (mna, N, subs_g)


def expand_matrix(mat, add_a_row=False, add_a_col=False):
    if add_a_row:
        row = sympy.zeros((1, mat.shape[1]))
        mat = mat.row_insert(mat.shape[0], row)
    if add_a_col:
        col = sympy.zeros((mat.shape[0], 1))
        mat = mat.col_insert(mat.shape[1], col)
    return mat


def get_roots(expr):
    num, den = sympy.fraction(expr)
    return sympy.solve(den, s), sympy.solve(num, s)


def parse_substitutions(slist):
    """Generates a substitution dictionary to be passed to solve()

    **Parameters:**

    slist : a list of strings
        The elements of the list should be according to the syntax
        ``'<part_id1>=<part_id2>'``, eg ``'R2=R1'``, instructing the simulator
         to use the value of R1 (R1) instead of R2.

    **Returns:**

    subs : dict
        the dictionary of symbols to be passed to :func:`solve`.
    """
    subs = {}
    for l in slist:
        v1, v2 = l.split("=")
        letter_id1 = v1[0].upper() if v1[0].upper() != 'R' else 'G'
        letter_id2 = v2[0].upper() if v2[0].upper() != 'R' else 'G'
        if letter_id1[0] in ('R', 'G', 'L', 'C', 'M'):
            s1 = _symbol_factory(letter_id1 + v1[1:], real=True, positive=True)
        else:
            s1 = _symbol_factory(letter_id1 + v1[1:], real=True)
        if letter_id2[0] in ('R', 'G', 'L', 'C', 'M'):
            s2 = _symbol_factory(letter_id2 + v2[1:], real=True, positive=True)
        else:
            s2 = _symbol_factory(letter_id2 + v2[1:], real=True)
        subs.update({s1:s2})
    return subs

def _symbol_factory(name, **options):
    filtered_options = {}
    for i in options:
        if options[i] and enabled_assumptions[i]:
            filtered_options.update({i:options[i]})
        else:
            pass # discarded
    return sympy.Symbol(name, **filtered_options)
