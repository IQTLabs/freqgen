from collections import defaultdict
from numpy.random import choice
from CAI import genetic_codes

def amino_acid_seq(length: int, frequencies: dict) -> str:
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
        sequence += choice(amino_acids, p=frequencies)
    return sequence

def amino_acids_to_codons(amino_acids: str, codon_frequencies: dict, genetic_code=1) -> str:
    '''Generates a DNA representation of an amino acid sequence.

    Args:
        amino_acids (str): The amino acids to convert to DNA.
        codon_frequencies (dict): A dictionary of codon frequencies for each amino acid. For each amino acid, the sum of the frequencies of its codon must be 1.
        genetic_code (int, optional): The genetic code to use when converting to DNA. Defaults to 1, the standard genetic code.

    Returns:
        str: A DNA sequence with the given codon usage.
    '''

    # create a translation table
    translation_table = defaultdict(list)
    for key, value in genetic_codes[genetic_code].items():
    	translation_table[value].append(key)
    translation_table = dict(translation_table)

    # generate the sequence
    sequence = ""
    for aa in amino_acids:
        codons_for_aa = translation_table[aa]
        sequence += choice(codons_for_aa, p=[codon_frequencies[codon] for codon in codons_for_aa])

    return sequence

def gc_content(seq):
    """Calculates the GC content of a sequence.

    Args:
        seq (str): The DNA sequence whose GC content is being calculated.

    Returns:
        float: The GC content.

    Raises:
        ValueError: When there is an invalid character in the sequence.
    """
    for i in seq:
        if i not in ["A", "T", "G", "C"]:
            raise ValueError("Invalid character in sequence.")
	return (seq.count("G") + seq.count("C")) / len(seq)
