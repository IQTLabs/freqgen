CLI Reference
=============

The output of the ``--help`` commands for the various Freqgen CLI tools is
included below for reference.

To summarize the relationship between the different commands, here's a graph of
the commands, the file types they take, and how they relate:

.. mermaid::

  graph TD
  aa[freqgen aa]-->|.faa|generate[freqgen generate]
  seq["DNA sequence(s)"]-->|.fna |featurize[freqgen featurize]
  featurize-->|.yaml|generate
  seq["DNA sequence(s)"]-->|.fna |aa
  amino["Amino acid sequence(s)"]-->|.faa|aa
  generate-->|.fna|vis[freqgen visualize]
  amino-->|.faa|generate

Overview
--------

.. command-output:: freqgen --help

Generation Reference
--------------------

.. command-output:: freqgen generate --help

Amino Acid Generation Reference
-------------------------------

.. command-output:: freqgen aa --help

Sequence Featurization Reference
--------------------------------

.. command-output:: freqgen featurize --help

Result Visualization Reference
------------------------------

.. command-output:: freqgen visualize --help
