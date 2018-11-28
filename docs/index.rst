.. freqgen documentation master file, created by
   sphinx-quickstart on Fri May  4 22:46:48 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Freqgen
=======

|Build Status| |CodeFactor| |Docs|

`Freqgen <https://github.com/Lab41/freqgen>`_ is a tool to generate coding DNA
sequences with specified amino acid usage frequencies or sequence, GC content,
codon usage bias, and/or :math:`k`-mer usage bias. To accomplish this, Freqgen
uses genetic algorithms to efficiently search the solution space of possible DNA
sequences to find ones that most closely match the desired parameters.

Features
--------

- Supports both CLI and Python module usage
- Thoroughly documented with examples
- Leverages NumPy for C-optimized number crunching
- Can simultaneously match multiple DNA statistics

Installation
------------

Simply run::

$ pip install freqgen

Or, to get the latest (but not necessarily stable) development version::

$ pip install git+https://github.com/Lab41/freqgen.git

Five-second CLI tutorial
------------------------

The basic flow of Freqgen can be summarized in three steps:

#. Generate a new amino acid sequence based on the amino acid usage profile of reference sequences. If you already have a specific amino acid sequence in mind (*i.e.* for synthetic biology uses), skip this step::

    $ freqgen aa reference_sequences.fna -o new_sequence.faa -l LENGTH

#. Create a YAML file containing :math:`k`-mer frequencies for the amino acid sequence's DNA to have::

    $ freqgen featurize reference_sequences.fna -k INT -o reference_freqs.yaml

#. Generate the DNA sequence coding for the amino acid sequence::

    $ freqgen -t reference_freqs.yaml -s new_sequence.faa -v -o optimized.fna

#. Visualize the results of the optimization (*optional*)::

    $ freqgen visualize --target reference_freqs.yaml --optimized optimized.fna


Citation
--------

To be determined.


.. |Build Status| image:: https://travis-ci.org/Lab41/freqgen.svg?branch=master
   :target: https://travis-ci.org/Lab41/freqgen

.. |Cov| image:: https://codecov.io/gh/Lab41/freqgen/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/Lab41/freqgen

.. |Docs| image:: http://readthedocs.org/projects/freqgen/badge/?version=latest
   :target: http://freqgen.readthedocs.io/en/latest/?badge=latest

.. |CodeFactor| image:: https://www.codefactor.io/repository/github/Lab41/freqgen/badge
   :target: https://www.codefactor.io/repository/github/Lab41/freqgen/


Table of Contents
-----------------

.. toctree::
   :maxdepth: 2

   methods
   usage
   api
   cli

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
