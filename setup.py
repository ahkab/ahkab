#!/usr/bin/env python
from distutils.core import setup
from ahkab.ahkab import __version__

setup(name='ahkab',
      version=__version__,
      description='Ahkab: a SPICE-like electronic circuit simulator',
      author='Giuseppe Venturini and others',
      author_email='ggventurini+GITHUB@gmail.com',
      url='http://ahkab.github.io/ahkab/',
      packages=['ahkab']
     )

print """
+---------------------------------------------------+
The following dependencies are needed to run ahkab:
- numpy: http://numpy.scipy.org/ 
- sympy: http://code.google.com/p/sympy/
- matplotlib: http://matplotlib.sourceforge.net/

They are available through PyPi. (See Install.md)
+---------------------------------------------------+
"""
