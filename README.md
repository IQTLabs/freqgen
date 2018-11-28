# Freqgen

<p align ="center">
<img src='https://raw.githubusercontent.com/Lab41/freqgen/master/logo/Freqgen2-01_icon_only.png' height="150">
</p>

[![Build Status](https://travis-ci.org/Lab41/freqgen.svg?branch=master)](https://travis-ci.org/Lab41/freqgen) [![Documentation Status](https://readthedocs.org/projects/freqgen/badge/?version=latest)](https://freqgen.readthedocs.io/en/latest/?badge=latest) [![CodeFactor](https://www.codefactor.io/repository/github/lab41/freqgen/badge)](https://www.codefactor.io/repository/github/lab41/freqgen)


Freqgen is a tool to generate coding DNA sequences with specified amino acid
usage frequencies or sequence, GC content, codon usage bias, and/or *k*-mer
usage bias. To accomplish this, Freqgen uses genetic algorithms to efficiently
search the solution space of possible DNA sequences to find ones that most
closely match the desired parameters.

## Features

- CLI and Python API
- Can simultaneously match multiple DNA statistics
- Built-in visualization utility
- Supports several fitness metrics (and you can bring your own!)

## Installation

Simply run:

    $ pip install freqgen

Or, to get the latest (but not necessarily stable) development version:

    $ pip install git+https://github.com/Lab41/freqgen.git

## Five-second CLI tutorial

The basic flow of Freqgen can be summarized in three steps:

1. Generate a new amino acid sequence based on the amino acid usage profile of
reference sequences. If you already have a specific amino acid sequence in mind
(*i.e.* for synthetic biology uses), skip this step:

        $ freqgen aa reference_sequences.fna -o new_sequence.faa -l LENGTH

2. Create a YAML file containing *k*-mer frequencies for the amino acid
sequence's DNA to have:

        $ freqgen featurize reference_sequences.fna -k INT -o reference_freqs.yaml

3. Generate the DNA sequence coding for the amino acid sequence:

        $ freqgen -t reference_freqs.yaml -s new_sequence.faa -v -o optimized.fna

4. Visualize the results of the optimization (*optional*):

        $ freqgen visualize --target reference_freqs.yaml --optimized optimized.fna

## Documentation

Read the full docs over at
[freqgen.readthedocs.io](http://freqgen.readthedocs.io).

## Citation

To be determined!
