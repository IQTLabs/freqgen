Methods
=======

This section covers the nitty-gritty details of the method by which Freqgen
operates. To get to the details of how to use the software to generate new DNA
sequences, feel free to skip ahead to the :ref:`Usage` documentation page.

Amino Acid Sequence Generation
------------------------------

To generate a novel amino acid sequence, we require a frequency for each amino
acid which can either be supplied directly or calculated from a set of reference
sequences. Then, an amino acid sequence of the desired length is constructed
from amino acids chosen with a probability proportional to its frequency.

To illustrate this method, suppose that the desired length is 12 amino acids and
a reference set of proteins is :math:`\{QA, QQ\}`. From this, we can see that
two amino acids, ``A`` and ``Q``, occur with 25% and 75% frequency in the
reference set, respectively. We then choose from :math:`\{A, Q\}` with a
probability of :math:`\{0.25, 0.75\}`. This process is repeated a total of 12
times. Thus, a potential new amino acid sequence would be ``AAQAQQQQQAQQ``. As
expected, ``A`` comprised 25% of the amino acids in the sequence with the other
75% being ``Q``. This method preserves the overall distribution of the original
distribution while having the advantage of being both fast and scaling linearly
with the desired length.

:math:`k`-mer Optimization
--------------------------

The method by which we use genetic algorithms to create DNA representations of
amino acid sequences that have a desired :math:`k`-mer usage bias can be
summarized as:

#. Given target sequence(s), calculate their :math:`k`-mer and/or codon usage biases.
#. Create a population of :math:`n` candidate synonymous sequences at random.
#. Measure each candidate sequence's fitness with respect to the target frequencies.
#. Create a new generation of synonymous solutions based on the fittest sequences in the last generation, carrying the best solution forward unchanged to prevent a decline in fitness between generations.
#. For each individual in the new generation, with probability :math:`m`, randomly "mutate" it while preserving synonymity.
#. Similarly, with probability :math:`c`, recombine two of the new sequences to produce two new synonymous sequences.
#. Repeat the process until the fitness of the fittest member of the population has stopped improving for :math:`g` generations.

where :math:`k`, :math:`n`, :math:`m`, :math:`c`, and :math:`g` are
user-configurable parameters.

Fitness Calculation
___________________

Freqgen supports two similar metrics of candidate sequence fitness: `Euclidean
distance <https://en.wikipedia.org/wiki/Euclidean_distance>`_ and
`Jensen-Shannon divergence
<https://en.wikipedia.org/wiki/Jensen–Shannon_divergence>`_ as implemented in
`dit <https://github.com/dit/dit>`_. By default, Freqgen uses Euclidean distance
for performance reasons. Jensen-Shannon divergence is provided as an alternative
fitness metric due the fact that it is bounded by :math:`[0,1]` whereas
Euclidean distance is unbounded and dependent on :math:`k`, therefore making
cross-optimization comparison difficult.

Mutation and Crossover
______________________

To mutate a candidate solution sequence, we choose a codon as well as a
different synonymous codon for it, both with equal probability. We then replace
the codon with its synonym and return the mutated synonymous sequence.

For crossover between two parents, we choose a codon index at random. Then, we
perform synonymous recombination by swapping the sequences beyond that index.
For example, if the two sequences for protein sequence ``MDT*`` were
``ATG-GAT-ACT-TGA`` and ``ATG-GAC-ACA-TAA`` and the zero-based index 1 were
chosen, then the two resultant sequences would be ``ATG-GAT-ACA-TAA`` and
``ATG-GAC-ACT-TGA``. Note that both of these sequences have the same amino acid
sequence.

This approach has a significant advantage over the naive approach of randomly
changing one nucleotide in the sequence that a traditional genetic algorithm
approach based on real world polymorphisms might use. By replacing a randomly
chosen codon with a synonymous one, we can ensure that the mutated individual
has the same amino acid sequence as before, eliminating the need to both verify
amino acid sequence identity and repeat the mutation function again until a
synonymous sequence is generated, as a naive nucleotide swap would require.
