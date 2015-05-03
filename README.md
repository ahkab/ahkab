<img src="https://rawgithub.com/ahkab/ahkab/master/doc/images/logo_small.png" alt="Monkeying around" style="width: 80px;"/> ahkab
=================================================================================================================================

**a SPICE-like electronic circuit simulator written in Python**

The code should be easy to read and modify, the main language is Python -- 2 or 3 -- and it is platform-independent.

## Welcome to the PYPY branch of `ahkab`!

* We try to run always against the latest [pypy](http://pypy.org/).
* Right now, most efforts are geared towards the Python 2 implementation, but that's just temporary.

All which is not strictly necessary is disabled in this experimental branch.

## Performance

Although numeric code does not really work for now with Pypy (we need `scipy`!), symbolic calculations are performed
very quickly on this alternative Python implementation.

Take a look at PYPY running 5+ times faster than CPython with no speed optimization on our side:


```
          Runtime [s]
          0.000     1.295     2.590     3.885     5.180     6.475     7.770     9.065     10.360     
          |         |         |         |         |         |         |         |         |         
 Nodes    +--------------------------------------------------------------------------------->
          |
    2.000 +o
          |=*
          |
    4.000 +o
          |==*
          |
    6.000 +o
          |===*
          |
    8.000 +----o
          |=====*
          |
   10.000 +--------------------o
          |=======*
          |
   12.000 +-------------------------------------------------------------------------------o
          |============*
          |
          +--------------------------------------------------------------------------------->
          |         |         |         |         |         |         |         |         |         
          0.000     1.295     2.590     3.885     5.180     6.475     7.770     9.065     10.360     
          Runtime [s]

    LEGEND:
    --* Ahkab 0.16-pypy, running on CPython, R2R ladder DAC
    ==o Ahkab 0.16-pypy, running on Pypy2, R2R ladder DAC
```
