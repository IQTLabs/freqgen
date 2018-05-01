from numpy.random import choice

def amino_acid_seq(length, frequencies) -> str:
    """Generates an amino acid sequence given frequencies of each amino acid.

    Args:
        length (int): The length of the amino acid sequence to generate.
        frequencies (dict): A dictionary containing a mapping of each amino acid to its frequency.

    Note:
        The sum of all the values in `frequencies` must be 1.

    Returns:
        str: An amino acid sequence with the given frequencies.
    """

    if length <= 0:
        raise ValueError("Length must be a positive integer.")

    if sum(frequencies.values()) != 1:
        raise ValueError("Sum of frequencies must be 1.")

    sequence = ""
    amino_acids, frequencies = zip(*frequencies.items())
    for i in range(length):
        sequence += choice(amino_acids, p=frequencies)
    return sequence
