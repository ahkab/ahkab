Installing ahkab
''''''''''''''''

Requirements
------------

The program requires:

-  the **Python 2** or **Python 3** interpreter (at least v.2.6 for
   Python2, at least v.3.3 for Python3),
-  **numpy>=1.7.0**,
-  **scipy>=0.14.0**,
-  **sympy>=0.7.5**,
-  and **tabulate>0.7.3**.

Strongly recommended:

-  **matplotlib>=1.1.1**,
-  **nose** for running the test suite.

Please try to use an up-to-date version of the libraries instead of the
bare minimum required.

All platform that are supported by the dependencies are also platforms
supported by ahkab, although the author only runs \*UNIX variants. If
you run into any problem, please report it in `the issue manager
<https://github.com/ahkab/ahkab/issues>`__.

--------------

`Numpy and Scipy <http://www.scipy.org/>`__ are needed for all the
numeric computations. On a Debian system, Python, Numpy and Scipy may be
installed running:

``# aptitude install python python-numpy python-scipy``

--------------

The symbolic analysis capabilities rely on the amazing
`sympy <http://www.sympy.org/>`__. Any version of sympy will do if you
are interested only in numeric simulations, but, if you run symbolic
simulations, *sympy version 0.7.5 or higher* is needed.

``# aptitude install python-sympy``

--------------

Plotting requires `matplotlib <http://matplotlib.sourceforge.net/>`__:

``# aptitude install python-matplotlib``

Install
-------

The project is hosted on GitHub and `on
PyPI <https://pypi.python.org/pypi/ahkab/>`__: to run ahkab, you can:

- run ``pip install ahkab --user``,
- `download a tarball containing the source
  code <https://github.com/ahkab/ahkab/archive/master.zip>`__,
- or check out the latest code `as explained here
  <https://help.github.com/articles/fetching-a-remote/#clone>`__.

If you choose not to use PIP, you will need to install the module
manually for your user on your \*UNIX system with the ``distutils``
script provided:

``python setup.py install --prefix=~/.local``

Thanks
------

Many thanks to the developers of the above libraries, their effort made
this project possible. :)
