from collections import defaultdict, Counter
from itertools import islice, product
import numpy as np
from CAI import genetic_codes
from warnings import warn

def amino_acid_seq(length, frequencies):
    """Generates an amino acid sequence given frequencies of each amino acid.

    Args:
        length (int): The length of the amino acid sequence to generate.
        frequencies (dict): A dictionary containing a mapping of each amino acid to its frequency.

    Note:
        The sum of all the values in ``frequencies`` must be 1.

    Returns:
        str: An amino acid sequence with the given frequencies.

    Raises:
        ValueError: When the length of the sequence is invalid or when the probabilities do not sum to 1.
    """

    if length <= 0:
        raise ValueError("Length must be a positive integer")

    sequence = ""
    amino_acids, frequencies = zip(*frequencies.items())
    for i in range(length):
        sequence += np.random.choice(amino_acids, p=frequencies)
    return sequence

def amino_acids_to_codons(aa_seq, codon_frequencies, genetic_code=1):
    '''Generates a DNA representation of an amino acid sequence.

    Args:
        aa_seq (str): The amino acids to convert to DNA.
        codon_frequencies (dict): A dictionary of codon frequencies for each amino acid. For each amino acid, the sum of the frequencies of its codon must be 1.
        genetic_code (int, optional): The genetic code to use when converting to DNA. Defaults to 1, the standard genetic code.

    Returns:
        str: A DNA sequence with the given codon usage.
    '''

    codons_dict = codons_for_aa(genetic_code)

    # generate the sequence
    sequence = ""
    for aa in aa_seq:
        try:
            codons = codons_dict[aa]
            sequence += np.random.choice(codons, p=[codon_frequencies[codon] for codon in codons])
        except KeyError:
            pass

    return sequence

def gc_content(dna_seq):
    """Calculates the GC content of a sequence.

    Args:
        dna_seq (str): The DNA sequence whose GC content is being calculated.

    Returns:
        float: The GC content.

    Raises:
        ValueError: When there is an invalid character in the sequence, i.e. not A, T, G, or C.
    """
    for i in dna_seq:
        if i not in ["A", "T", "G", "C"]:
            raise ValueError("Invalid character in sequence.")
    return (dna_seq.count("G") + dna_seq.count("C")) / len(dna_seq)

def codons_for_aa(genetic_code):
    '''Generates a dict of the codons for each amino acid.

    Args:
        genetic_code (int): The genetic code to use to create the dictionary.

    Returns:
        dict: A dict with the amino acids as keys and a list of codons for the amino acid as the values.

    Example:
        >>> codons_for_aa(1)
        {'A': ['GCT', 'GCC', 'GCA', 'GCG'],
         'C': ['TGT', 'TGC'],
         'D': ['GAT', 'GAC'],
         'E': ['GAA', 'GAG'],
         'F': ['TTT', 'TTC'],
         'G': ['GGT', 'GGC', 'GGA', 'GGG'],
         'H': ['CAT', 'CAC'],
         'I': ['ATT', 'ATC', 'ATA'],
         'K': ['AAA', 'AAG'],
         'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
         'M': ['ATG'],
         'N': ['AAT', 'AAC'],
         'P': ['CCT', 'CCC', 'CCA', 'CCG'],
         'Q': ['CAA', 'CAG'],
         'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
         'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
         'T': ['ACT', 'ACC', 'ACA', 'ACG'],
         'V': ['GTT', 'GTC', 'GTA', 'GTG'],
         'W': ['TGG'],
         'Y': ['TAT', 'TAC']}
    '''
    # create a translation table
    codons_for_aa = defaultdict(list)
    for key, value in genetic_codes[genetic_code].items():
    	codons_for_aa[value].append(key)
    return dict(codons_for_aa)

def codon_frequencies(dna_seq, genetic_code=1):
    '''Calculated the codon frequencies of each codon within the amino acid it represents.

    Args:
        dna_seq (str): The DNA sequence.
        genetic_code (int, optional): The genetic code to use. Defaults to the standard genetic code.

    Returns:
        dict: The codon frequencies of each codon.

    Raises:
        ValueError: When the sequence length is not divisible into codons, i.e. when sequence length is not divisible by three.
    '''

    if len(dna_seq) % 3 != 0:
        raise ValueError("Sequence length is not divisible by 3.")

    dna_seq = [dna_seq[i:i+3] for i in range(0, len(dna_seq), 3)]
    counts = Counter(dna_seq)

    freqs = dict()
    for aa, codons in codons_for_aa(genetic_code).items():
        total_codons_for_aa = sum([counts[codon] for codon in codons])
        for codon in codons:
            try:
                freqs[codon] = counts[codon] / total_codons_for_aa
            except ZeroDivisionError:
                freqs[codon] = 0
    return freqs

def translate(dna_seq, genetic_code=1):
    """Translates a DNA sequence into amino acids.

    Args:
        dna_seq (str): The DNA sequence to translate.
        genetic_code (int, optional): The genetic code to use. Defaults to the standard genetic code.

    Note:
        If there is a stop codon in the sequence, it is ignored.

    Returns:
        str: The amino acid sequence.

    Raises:
        ValueError: When there is an invalid character or when the sequence
            length is not divisible into codons, i.e. when sequence length is not
            divisible by three.

    Example:
        >>> translate("ATTAATCAAACGGAGTTA")
        'INQTEL'
    """

    for base in dna_seq:
        if base not in ["A", "T", "G", "C"]:
            raise ValueError("Invalid character in sequence: ", base)
    if len(dna_seq) % 3 != 0:
        raise ValueError("Invalid sequence length.")

    codons = [dna_seq[i:i+3] for i in range(0, len(dna_seq), 3)]
    aa_seq = ""
    for i, codon in enumerate(codons):
        try:
            aa_seq += genetic_codes[genetic_code][codon]
        except KeyError:
            warn("Stop codon in sequence... Ignoring!")
            pass
    return aa_seq

def k_mers(seq, k):
    '''Yields all *k*-mers in the input sequence with repeats.

    Modified from this `StackOverflow answer <https://stackoverflow.com/questions/6822725/rolling-or-sliding-window-iterator>`_.

    Args:
        seq (str): The sequence for which to generate *k*-mers.
        k (int): the length of the *k*-mers.

    Yields:
        str: the next *k*-mer

    Example:
        >>> list(k_mers("GATTACA", 1))
        ['G', 'A', 'T', 'T', 'A', 'C', 'A']
        >>> list(k_mers("GATTACA", 2))
        ['GA', 'AT', 'TT', 'TA', 'AC', 'CA']
        >>> list(k_mers("GATTACA", 3))
        ['GAT', 'ATT', 'TTA', 'TAC', 'ACA']
        >>> list(k_mers("GATTACA", 4))
        ['GATT', 'ATTA', 'TTAC', 'TACA']
        >>> k_mers("GATTACA", 4)
        <generator object k_mers at 0x10831d258>
    '''
    it = iter(seq)
    result = tuple(islice(it, k))
    if len(result) == k:
        yield "".join(result)
    for elem in it:
        result = result[1:] + (elem,)
        yield "".join(result)

def k_mer_frequencies(seq, k):
    '''Calculates relative frequencies of each *k*-mer in the sequence.

    Args:
        seq (str): The sequence to for which to generate *k*-mer frequencies.
        k (int): the length of the *k*-mers.

    Returns:
        dict: A dict in which the keys are *k*-mers and the values are floats of their frequencies.

    Example:
        >>> k_mer_frequencies("INQTEL", 1)
        {'E': 0.16666666666666666,
         'I': 0.16666666666666666,
         'L': 0.16666666666666666,
         'N': 0.16666666666666666,
         'Q': 0.16666666666666666,
         'T': 0.16666666666666666}

         >>> k_mer_frequencies("GATGATGGC", 3)
         {'ATG': 0.2857142857142857,
          'GAT': 0.2857142857142857,
          'GGC': 0.14285714285714285,
          'TGA': 0.14285714285714285,
          'TGG': 0.14285714285714285}
    '''

    count = Counter(k_mers(seq, k))
    return {k: v/sum(count.values()) for k, v in count.items()}
