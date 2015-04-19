#!/usr/bin/env python
import os
import re
from setuptools import setup, find_packages

with open('ahkab/__version__.py') as version_file:
    __version__ = re.findall("__version__ = '([^']*)'", version_file.read())[0]

def read(fname):
    try:
        with open(os.path.join(os.path.dirname(__file__), fname)) as fp:
            return fp.read()
    except IOError:
        return ""

setup(
    name='ahkab',
    version=__version__,
    packages=find_packages(),
    package_data={'ahkab': ['test_data/*.dat']},
    install_requires=['numpy>1.7.0', 'scipy>=0.14.0', 'sympy>=0.7.6',
                      'tabulate>=0.7.3'],
    extras_require = {'PLOT':['matplotlib>=1.1.1']},
    zip_safe=False,
    include_package_data=True,
    author="Giuseppe Venturini and others",
    author_email="ggventurini+GITHUB@gmail.com",
    description="a SPICE-like electronic circuit simulator",
    long_description=''.join([read('pypi_description.rst'), '\n\n',
                              read('CHANGES.rst')]),
    entry_points = {'console_scripts': ['ahkab = ahkab.__main__:_cli', ], },
    license="GPL",
    keywords="electronic circuit simulator numeric symbolic",
    url='http://ahkab.github.io/ahkab/',
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Operating System :: POSIX",
        "Operating System :: POSIX :: Linux",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: MacOS",
        "Natural Language :: English",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.4"]
)

