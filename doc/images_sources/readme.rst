This directory contains the TeX sources for the circuit symbols
used in the documentation.

To generate the SVG files, you need a latex distribution, tikz/pgf, circuitikz,
and the utility pdf2svg.

Eg.

::

    pdflatex resistor.tex
    pdf2svg resistor.pdf resistor.svg

Many thanks to the authors of the above packages and tools.
