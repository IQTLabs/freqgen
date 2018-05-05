from collections import defaultdict, Counter
import numpy.random
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
        ValueError: When the length of the sequence is invalid.
        ValueError: When the probabilities do not sum to 1.
    """

    if length <= 0:
        raise ValueError("Length must be a positive integer.")

    sequence = ""
    amino_acids, frequencies = zip(*frequencies.items())
    for i in range(length):
        sequence += numpy.random.choice(amino_acids, p=frequencies)
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
            sequence += numpy.random.choice(codons, p=[codon_frequencies[codon] for codon in codons])
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
        ValueError: When there is an invalid character in the sequence.
    """
    for i in dna_seq:
        if i not in ["A", "T", "G", "C"]:
            raise ValueError("Invalid character in sequence.")
    return (dna_seq.count("G") + dna_seq.count("C")) / len(dna_seq)

def codons_for_aa(genetic_code):
    '''Generates a dict of the codons for each amino acid.

    Example:
        >>> codons_for_aa(1)
        {'F': ['TTT', 'TTC'], 'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'Y': ['TAT', 'TAC'], 'C': ['TGT', 'TGC'], 'W': ['TGG'], 'P': ['CCT', 'CCC', 'CCA', 'CCG'], 'H': ['CAT', 'CAC'], 'Q': ['CAA', 'CAG'], 'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'I': ['ATT', 'ATC', 'ATA'], 'M': ['ATG'], 'T': ['ACT', 'ACC', 'ACA', 'ACG'], 'N': ['AAT', 'AAC'], 'K': ['AAA', 'AAG'], 'V': ['GTT', 'GTC', 'GTA', 'GTG'], 'A': ['GCT', 'GCC', 'GCA', 'GCG'], 'D': ['GAT', 'GAC'], 'E': ['GAA', 'GAG'], 'G': ['GGT', 'GGC', 'GGA', 'GGG']}
    '''
    # create a translation table
    codons_for_aa = defaultdict(list)
    for key, value in genetic_codes[genetic_code].items():
    	codons_for_aa[value].append(key)
    return dict(codons_for_aa)

def codon_frequencies(dna_seq, genetic_code=1):
    '''Calculated the codon frequencies of each codon

    Args:
        dna_seq (str): The DNA sequence.
        genetic_code (int, optional): The genetic code to use. Defaults to the standard genetic code.

    Returns:
        dict: The codon frequencies of each codon.
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
    """Translates a DNA sequence.
    """
    codons = [dna_seq[i:i+3] for i in range(0, len(dna_seq), 3)]
    aa_seq = ""
    for i, codon in enumerate(codons):
        try:
            aa_seq += genetic_codes[genetic_code][codon]
        except KeyError:
            warn("Stop codon in sequence... Ignoring!")
            pass
    return aa_seq
