CLI Reference
=============

The output of the ``--help`` commands for the various Freqgen CLI tools is
included below for reference.

Overview
--------

.. code::

    $ freqgen --help
    Usage: freqgen [OPTIONS] COMMAND [ARGS]...

    Options:
      --help  Show this message and exit.

    Commands:
      generate*  Generate a new DNA sequence with matching...
      aa         Generate an amino acid sequence from FASTA
      featurize  Featurize a FASTA file

Generation Reference
--------------------

.. code::

    $ freqgen generate --help
    Usage: freqgen generate [OPTIONS]

      Generate a new DNA sequence with matching features

    Options:
      -s, --seq PATH             The target amino acid sequence.
      -f, --freqs PATH           The target frequencies.
      -v, --verbose              Whether to show optimization progress. Defaults
                                 to false.
      -i INTEGER                 How many generations to stop after no
                                 improvement. Defaults to 50.
      -p INTEGER                 Population size. Defaults to 100.
      -m FLOAT                   Mutation rate. Defaults to 0.3.
      -c FLOAT                   Crossover rate. Defaults to 0.8.
      -t, --trans-table INTEGER  The translation table to use. Defaults to 11, the
                                 standard genetic code.
      -o, --output PATH
      --help                     Show this message and exit.

Amino Acid Generation Reference
-------------------------------

.. code::

    $ freqgen aa --help
    Usage: freqgen aa [OPTIONS] FILEPATH

      Generate an amino acid sequence from FASTA

    Options:
      --mode [freq|seq]          Whether to use the exact AA seq or its
                                 frequencies. Defaults to freq.
      -t, --trans-table INTEGER  The translation table to use. Defaults to 11, the
                                 standard genetic code.
      -l, --length INTEGER       The length of the AA sequence (excluding stop
                                 codon) to generate if --mode=freq.
      -s, --stop-codon           Whether to include a stop codon. Defaults to
                                 true.
      -v, --verbose              Whether to print final result if outputting to
                                 file. Defaults to false.
      -o, --output PATH          The output file.
      --help                     Show this message and exit.

Sequence Featurization Reference
--------------------------------

.. code::

    $ freqgen featurize --help
    Usage: freqgen featurize [OPTIONS] FILEPATH

      Featurize a FASTA file

    Options:
      -k INTEGER         Values of k to featurize the seqs for. May be repeated.
      -c, --codon-usage  Whether to include a codon frequency featurization.
      -o, --output PATH  The output file.
      --help             Show this message and exit.
