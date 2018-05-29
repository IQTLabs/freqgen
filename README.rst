Freqgen
=======
|Build Status| |Docs|

Freqgen is a tool to generate coding DNA sequences with specified GC content,
amino acid usage, codon usage bias, and *k*-mer usage bias. Freqgen uses genetic
algorithms to efficiently search the solution space of possible DNA sequences to
find sequences that most closely match the desired parameters.

Features
--------

- Supports both CLI and Python module usage
- Thoroughly documented with examples
- Leverages NumPy for C optimized number crunching
- Can simultaneously match multiple DNA statistics

Installation
------------

Simply use::

$ pip install freqgen

Or, to get the newest development version::

$ git clone https://github.com/benjamin-lee/freqgen.git
$ cd freqgen
$ python setup.py install

Citation
--------

To be determined.


.. |Build Status| image:: https://travis-ci.org/Benjamin-Lee/freqgen.svg?branch=master
   :target: https://travis-ci.org/Benjamin-Lee/freqgen

.. |Docs| image:: http://readthedocs.org/projects/freqgen/badge/?version=latest
   :target: http://freqgen.readthedocs.io/en/latest/?badge=latest
