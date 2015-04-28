Installing ahkab
''''''''''''''''

Requirements
------------

The program requires:

-  the **Python 2** or **Python 3** interpreter (at least v.2.6 for
   Python2, at least v.3.3 for Python3),
-  **numpy>=1.7.0**,
-  **scipy>=0.14.0**,
-  **sympy>=0.7.6**,
-  and **tabulate>0.7.3**.

Strongly recommended:

-  **matplotlib>=1.1.1**,
-  **nose** for running the test suite.

Please try to use an up-to-date version of the libraries instead of the
bare minimum required.

All platform that are supported by the dependencies are also platforms
supported by ``ahkab``, although the author only runs \*UNIX variants. If
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
simulations, *sympy version 0.7.6 or higher* is needed.

``# aptitude install python-sympy``

--------------

Plotting requires `matplotlib <http://matplotlib.sourceforge.net/>`__:

``# aptitude install python-matplotlib``

Install
-------

The source code for the project is `hosted on GitHub
<https://github.com/ahkab/ahkab>`__ and releases can be found `on PyPI
<https://pypi.python.org/pypi/ahkab/>`__.

To install ``ahkab``, you can have two options: using ``pip`` or using
``distutils``.

Install with pip
""""""""""""""""

If you use ``pip``, which boundled in your Python installation, the source code
is downloaded from you off the Python Package Index (`PyPI
<https://pypi.python.org/>`__).

You may:

- Issue ``pip install ahkab``, which may require administrative access depending
  on what permissions your user has on your Python installation.
- To avoid having to supply admin credentials, you may use ``pip`` according to
  "the user scheme", issuing ``pip install ahkab --user``.

Install with distutils
""""""""""""""""""""""

Installing manually through ``distutils`` requires that you download the source
code, untar and move to the root directory of the package.

For which you should first either:
  - `download a release tarball containing the source code
    <https://github.com/ahkab/ahkab/releases/>`__,
  - or check out the latest code `as explained on GitHub
    <https://help.github.com/articles/fetching-a-remote/#clone>`__.

Then, you will need to install the module manually with the ``distutils`` script ``setup.py`` provided, you can choose whether:
    - to install for all users: ``python setup.py install``
    - or only your own: ``python setup.py install --user``
    - or to install to a different prefix: ``python setup.py install
      --prefix=~/.local``

The Python documentation for installing with ``distutils`` will clear up any
remaining doubt: `for version 2 of the language
<https://docs.python.org/2/install/#the-new-standard-distutils>`__, `for version
3 <https://docs.python.org/3/install/>`__.

Thanks
------

Many thanks to the developers of the above libraries, their effort made
this project possible. :)
