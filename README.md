<img src="http://raw.github.com/wiki/ahkab/ahkab/images/logo_small.png" alt="Monkeying around" style="width: 80px;"/> ahkab
===========================================================================================================================

**a SPICE-like electronic circuit simulator written in Python**

The code should be easy to read and modify, the main language is Python -- 2 or 3 -- and it is platform-independent.

## Welcome to the PYPY branch of `ahkab`!

* We try to run always against the latest [pypy](http://pypy.org/).
* Right now, most efforts are geared towards the Python 2 implementation, but that's just temporary.

All which is not strictly necessary is disabled in this experimental branch.

## Performance

Although numeric code does not really work for now with Pypy (we need `scipy`!), symbolic calculations are performed
very quickly on this alternative Python implementation.

Take a look at Pypy running 5x faster than CPython with no speed optimization on our side:

```
      Runtime [s]
      0.000     1.505     3.009     4.514     6.018     7.523     9.028     10.532     12.037     
      |         |         |         |         |         |         |         |         |         
      +--------------------------------------------------------------------------------->
      |
2.000 +-o
      |*
      |
4.000 +---o
      |*
      |
6.000 +------o
      |=*
      |
8.000 +-------o
      |====*
      |
10.000 +-------o
      |====================*
      |
12.000 +----------------o
      |===============================================================================*
      |
      +--------------------------------------------------------------------------------->
      |         |         |         |         |         |         |         |         |         
      0.000     1.505     3.009     4.514     6.018     7.523     9.028     10.532     12.037
      Runtime [s]


LEGEND:
--o Ahkab 0.15-pypy, running on Pypy2, R2R ladder DAC
==* Ahkab 0.15-pypy, running on CPython, R2R ladder DAC
```
