# -*- coding: iso-8859-1 -*-
# test_fourier.py
# Unit tests for the FOURIER module
# Copyright 2015 Giuseppe Venturini
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
# along with ahkab. If not, see <http://www.gnu.org/licenses/>.

from __future__ import unicode_literals, print_function, division
import numpy as np
import ahkab

cir = ahkab.Circuit('Test FOUR and FFT')
#mys = lambda t: 1 if not t else math.sin(math.pi*1e4*t)/(math.pi*1e4*t)
mys = ahkab.time_functions.sin(vo=0, va=1, freq=10e3)
cir.add_resistor('R1', 'n1', cir.gnd, 1e3)
cir.add_vsource('V1', 'n1', cir.gnd, 1, function=mys)
tr = ahkab.new_tran(0, 1e-3, 1e-5, x0=None)
r = ahkab.run(cir, tr)['tran']



def test_fourier():
    """Test fourier.fourier() and printing.print_fourier()"""
    fs, F, THD = ahkab.fourier.fourier('vn1', r, 10e3)
    fs_ref = np.array([0., 10000., 20000., 30000., 40000., 50000., 60000.,
                       70000., 80000., 90000.])
    assert np.allclose(fs, fs_ref)
    assert np.argmax(abs(F)) == 1
    assert len(F) == len(fs) == 10
    assert abs(F[0]) < 1
    assert THD < .1
    ahkab.printing.print_fourier('vn1', fs, F, THD)

def test_spicefft():
    """Test fourier.spicefft() and printing.print_spicefft()"""
    for w in (ahkab.options.RECT_WINDOW, ahkab.options.BART_WINDOW,
              ahkab.options.HANN_WINDOW, ahkab.options.HAMM_WINDOW,
                                    ahkab.options.BLACK_WINDOW, ahkab.options.HARRIS_WINDOW,
                                    ahkab.options.KAISER_WINDOW):
        fsp, Fsp, THDsp = ahkab.fourier.spicefft('vn1', r, window=w)
    fsp, Fsp, THDsp = ahkab.fourier.spicefft('vn1', r, window='GAUSS', alpha=3000)
    ahkab.printing.print_spicefft('vn1', fsp, Fsp, uformat='UNORM', window='GAUSS')

    assert np.allclose(fsp, np.linspace(0, 512*1000, 512, endpoint=False))

    assert len(fsp) == 512
    assert len(Fsp) == 512

    fsp, Fsp, THDsp = ahkab.fourier.spicefft('vn1', r, window='GAUSS', alpha=3000, np=4096)
    assert len(fsp) == 2048
    assert len(Fsp) == 2048
    assert fsp[0] == 0
    assert np.real(Fsp[0]) == np.real_if_close(Fsp[0])

    f, F, THD = ahkab.fourier.spicefft('vn1', r, 10e3)
    ahkab.printing.print_spicefft('vn1', f, F, THD, outfile='stdout')

    f, F, THD = ahkab.fourier.spicefft('vn1', r, 10e3, fmin=20e3, fmax=50e3)
    assert np.allclose(f, np.array([0., 20000., 30000., 40000., 50000.]))
    ahkab.printing.print_spicefft('vn1', f, F, THD, outfile='stdout')

    f1, F1, THD1 = ahkab.fourier.spicefft('vn1', r, 10e3, window='HANN', **{'from':.1e-3, 'to':.4e-3})
    f2, F2, THD1 = ahkab.fourier.spicefft('vn1', r, 10e3, window='HANN', **{'start':.1e-3, 'stop':.4e-3})

    assert np.allclose(f1, f2)
    assert np.allclose(F1, F2)



