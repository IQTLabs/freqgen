from collections import defaultdict, Counter
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

    codons_dict = codons_for_aa(genetic_code)

    # generate the sequence
    sequence = ""
    for aa in amino_acids:
        codons = codons_dict[aa]
        sequence += choice(codons, p=[codon_frequencies[codon] for codon in codons])

    return sequence

def gc_content(seq: str) -> float:
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

def codons_for_aa(genetic_code: int) -> dict:
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

def codon_frequencies(seq: str, genetic_code=1) -> dict:
    '''Calculated the codon frequencies of each codon

    Args:
        seq (str): The DNA sequence.
        genetic_code (int, optional): The genetic code to use. Defaults to the standard genetic code.

    Returns:
        dict: The codon frequencies of each codon.
    '''

    if len(seq) % 3 != 0:
        raise ValueError("Sequence length is not divisible by 3.")

    seq = [seq[i:i+3] for i in range(0, len(seq), 3)]
    counts = Counter(seq)

    freqs = dict()
    for aa, codons in codons_for_aa(genetic_code).items():
        total_codons_for_aa = sum([counts[codon] for codon in codons])
        for codon in codons:
            try:
                freqs[codon] = counts[codon] / total_codons_for_aa
            except ZeroDivisionError:
                freqs[codon] = 0
    return freqs

def translate(seq, genetic_code=1):
    codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    aa_seq = ""
    for codon in codons:
        print(genetic_codes[genetic_code][codon])
    return aa_seq
with open("dummy_seq.txt") as f:
    # print(codon_frequencies(f.readline().rstrip()))
    print(translate(f.readline().rstrip()))
# print(codons_for_aa(1))
