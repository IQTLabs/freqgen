from collections import defaultdict, Counter, Iterable
from itertools import islice, product, chain
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

    Example:
        >>> from Bio import SeqIO
        >>> seq = SeqIO.read("beta_lactamase.fasta", "fasta").seq
        >>> frequencies = k_mer_frequencies(seq, 1)
        >>> amino_acid_seq(25, frequencies)
    """

    if length <= 0:
        raise ValueError("Length must be a positive integer")

    sequence = ""
    amino_acids, frequencies = zip(*frequencies.items())
    for i in range(length):
        sequence += np.random.choice(amino_acids, p=frequencies)
    return sequence

def amino_acids_to_codons(aa_seq, codon_frequencies, genetic_code=11):
    '''Generates a DNA representation of an amino acid sequence.

    Args:
        aa_seq (str): The amino acids to convert to DNA.
        codon_frequencies (dict): A dictionary of codon frequencies for each amino acid. For each amino acid, the sum of the frequencies of its codons must be 1.
        genetic_code (int, optional): The genetic code to use when converting to DNA. Defaults to 11, the standard genetic code.

    Returns:
        str: A DNA sequence with the given codon usage.

    Example:
        >>> from Bio import SeqIO
        >>> seq = SeqIO.read("sequence.fasta", "fasta").seq
        >>> frequencies = codon_frequencies(seq)
        >>> amino_acids_to_codons("INQTEL", frequencies)
        'ATAAATCAAACCGAACTT'
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

    Example:
        >>> gc_content("GATTACA")
        0.2857142857142857
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

def codon_frequencies(dna_seq, genetic_code=11):
    '''Calculates the codon frequencies of each codon within the amino acid it represents.

    Args:
        dna_seq (str): The DNA sequence.
        genetic_code (int, optional): The genetic code to use. Defaults to the 11, standard genetic code.

    Returns:
        dict: The codon frequencies of each codon.

    Raises:
        ValueError: When the sequence length is not divisible into codons, i.e. when sequence length is not divisible by three.

    Example:
        >>> from Bio import SeqIO
        >>> seq = SeqIO.read("sequence.fasta", "fasta").seq
        >>> codon_frequencies(seq)
        {'AAA': 0.5,
         'AAC': 0.5,
         'AAG': 0.5,
         'AAT': 0.5,
         'ACA': 0.0,
         'ACC': 0.5,
         'ACG': 0.0,
         'ACT': 0.5,
         'AGA': 0.0,
         'AGC': 0.0,
         'AGG': 0.0,
         'AGT': 0.0,
         'ATA': 1.0,
         'ATC': 0.0,
         'ATG': 1.0,
         'ATT': 0.0,
         'CAA': 1.0,
         'CAC': 0.0,
         'CAG': 0.0,
         'CAT': 1.0,
         'CCA': 0.25,
         'CCC': 0.125,
         'CCG': 0.375,
         'CCT': 0.25,
         'CGA': 0.3333333333333333,
         'CGC': 0.0,
         'CGG': 0.6666666666666666,
         'CGT': 0.0,
         'CTA': 0.0,
         'CTC': 0.0,
         'CTG': 0.25,
         'CTT': 0.75,
         'GAA': 1.0,
         'GAC': 0.0,
         'GAG': 0.0,
         'GAT': 1.0,
         'GCA': 0.0,
         'GCC': 0.8,
         'GCG': 0.2,
         'GCT': 0.0,
         'GGA': 0.5,
         'GGC': 0.25,
         'GGG': 0.25,
         'GGT': 0.0,
         'GTA': 0.0,
         'GTC': 0.25,
         'GTG': 0.0,
         'GTT': 0.75,
         'TAC': 0.0,
         'TAT': 1.0,
         'TCA': 0.2,
         'TCC': 0.4,
         'TCG': 0.0,
         'TCT': 0.4,
         'TGC': 1.0,
         'TGG': 1.0,
         'TGT': 0.0,
         'TTA': 0.0,
         'TTC': 0.6,
         'TTG': 0.0,
         'TTT': 0.4}
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

def translate(dna_seq, genetic_code=11):
    """Translates a DNA sequence into amino acids.

    Args:
        dna_seq (str): The DNA sequence to translate.
        genetic_code (int, optional): The genetic code to use. Defaults to 11, the standard genetic code.

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

    Raises:
        ValueError: When the value of *k* is less than the length of the sequence, k <= 0, or len(seq) is 0.

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

    # error checking
    if k > len(seq):
        raise ValueError("k may not be less then length of seq.")
    elif len(seq) == 0:
        raise ValueError("seq length may not be zero")
    elif k <= 0:
        raise ValueError("k may not be <= zero")

    it = iter(seq)
    result = tuple(islice(it, k))
    if len(result) == k:
        yield "".join(result)
    for elem in it:
        result = result[1:] + (elem,)
        yield "".join(result)

def k_mer_frequencies(seq, k, include_missing=False, vector=False):
    '''Calculates relative frequencies of each *k*-mer in the sequence.

    Args:
        seq (str or list): The sequence(s) to for which to generate *k*-mer frequencies.
        k (int or list): the length of the *k*-mer(s).
        include_missing (bool, optional): If True, include missing *k*-mers as having a frequency of 0. Only supports DNA *k*-mers. Defaults to False.
        vector (bool, optional): Return a 1-D Numpy array of the *k*-mer frequencies, ordered by *k*-mers alphabetically. If True, ``include_missing`` must also be True. Defaults to False.

    Returns:
        dict: A dict in which the keys are *k*-mers and the values are floats of their frequencies.

    Raises:
        ValueError: When an invalid value of k is provided or ``include_missing`` is False and ``vector`` is True.

    Example:
        >>> k_mer_frequencies("INQTEL", 1)
        {'E': 0.16666666666666666,
         'I': 0.16666666666666666,
         'L': 0.16666666666666666,
         'N': 0.16666666666666666,
         'Q': 0.16666666666666666,
         'T': 0.16666666666666666}

        >>> k_mer_frequencies("GATGATGGC", [1, 2])
        {'A': 0.2222222222222222,
         'AT': 0.25,
         'C': 0.1111111111111111,
         'G': 0.4444444444444444,
         'GA': 0.25,
         'GC': 0.125,
         'GG': 0.125,
         'T': 0.2222222222222222,
         'TG': 0.25}

        >>> k_mer_frequencies("GATGATGGC", 2, include_missing=True)
        {'AA': 0,
         'AC': 0,
         'AG': 0,
         'AT': 0.25,
         'CA': 0,
         'CC': 0,
         'CG': 0,
         'CT': 0,
         'GA': 0.25,
         'GC': 0.125,
         'GG': 0.125,
         'GT': 0,
         'TA': 0,
         'TC': 0,
         'TG': 0.25,
         'TT': 0}

        >>> k_mer_frequencies("GATGATGGC", 2, include_missing=True, vector=True)
        array([0.   , 0.   , 0.   , 0.25 , 0.   , 0.   , 0.   , 0.   , 0.25 ,
               0.125, 0.125, 0.   , 0.   , 0.   , 0.25 , 0.   ])
    '''

    if not include_missing and vector:
        raise ValueError("May not create vector without including missing kmers.")

    if not isinstance(k, Iterable):
        k = [k]
    else:
        k = sorted(k)

    output = []

    if isinstance(seq, str):
        seq = [seq]

    for _k in k:

        # check the value of k
        if _k < 1:
            raise ValueError("Invalid value of k. May not be less than 1.")

        # get all the k-mers for the seqs
        _seqs = []
        for _seq in [list(k_mers(_seq, _k)) for _seq in seq]:
            _seqs.extend(_seq)

        # determine their frequencies
        count = Counter(_seqs)
        frequencies = {key: value/sum(count.values()) for key, value in count.items()}

        if include_missing:
            defaults = {"".join(x): 0 for x in list(product("ATGC", repeat=_k))}
            frequencies = {**defaults, **frequencies}
        if vector:
            frequencies = sorted(list(frequencies.items()), key=lambda x: x[0])
            frequencies = np.fromiter((x[1] for x in frequencies), float, count=len(frequencies))
        output.append(frequencies)

    if len(output) == 1:
        return output[0]
    elif vector:
        return np.array(list(chain.from_iterable(output)))
    else:
        return {k: v for d in output for k, v in d.items() }

    return frequencies

