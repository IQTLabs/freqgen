import random
from warnings import warn
from math import isclose

import dit
import numpy as np
import yaml
from Bio.Seq import Seq
from dit.divergences import jensen_shannon_divergence

from .freqgen import *
from .pyeasyga import GeneticAlgorithm


def dna_to_vector(seq):
    seq = np.array(list(seq))
    seq[seq == "T"] = 0
    seq[seq == "C"] = 1
    seq[seq == "A"] = 2
    seq[seq == "G"] = 3
    seq = seq.astype(int)
    seq = [np.binary_repr(x, width=2) for x in seq]
    seq = [list(x) for x in seq]
    seq = np.array(seq).flatten().astype(int)
    return seq

def vector_to_dna(vector):
    dna = "".join(np.array(vector).astype(str))
    dna = np.array([dna[i:i + 2] for i in range(0, len(dna), 2)])
    dna[dna == "00"] = "T"
    dna[dna == "01"] = "C"
    dna[dna == "10"] = "A"
    dna[dna == "11"] = "G"
    return "".join(dna)

def _synonymous_codons(genetic_code_dict): # from CAI source code
    # invert the genetic code dictionary to map each amino acid to its codons
    codons_for_amino_acid = {}
    for codon, amino_acid in genetic_code_dict.items():
        codons_for_amino_acid[amino_acid] = codons_for_amino_acid.get(amino_acid, [])
        codons_for_amino_acid[amino_acid].append(codon)

    # create dictionary of synonmyous codons
    # Example: {'CTT': ['CTT', 'CTG', 'CTA', 'CTC', 'TTA', 'TTG'], 'ATG': ['ATG']...}
    return {codon: codons_for_amino_acid[genetic_code_dict[codon]] for codon in genetic_code_dict.keys()}

def generate(target_params,
             aa_seq,
             population_size=100,
             mutation_probability=0.3,
             crossover_probability=0.8,
             max_gens_since_improvement=50,
             improvement_rel_threshold=0.0,
             genetic_code=11,
             verbose=False,
             mode="ED",
             fitness_function=None):
    '''Generate a sequence matching :math:`k`-mer usage.

    Args:
        target_params (dict): The parameters to optimize towards. Should be of the format {:math:`k_n`: {:math:`k_{n1}`: 0.2, :math:`k_{n2}`: 0.3,...}...}. Pass absolute codon usage with ``"codons"`` as the key.
        aa_seq (str): The amino acid sequence for the optimized sequence.
        population_size (int, optional): The size of the population for the genetic algorithm. Defaults to 100.
        mutation_probability (float, optional): The likelihood of changing each member of each generation. Defaults to 0.3.
        crossover_probability (float, optional): The likelihood of each member of the population undergoing crossover. Defaults to 0.8.
        max_gens_since_improvement (int, optional): The number of generations of no improvement after which to stop optimization. Defaults to 50.
        improvement_rel_threshold (float, optional): The minimum ftness improvement in precentage for which to reset the value of generations since improvement. Defaults to 0, meaning that any improvement resets the counter. Greater values will result in earlier stopping.
        genetic_code (int, optional): The genetic code to use. Defaults to 11, the standard genetic code.
        verbose (bool, optional): Whether to print the generation number, generations since improvement, and fitness. Defaults to false.
        mode (str, optional): Whether to use Jensen-Shannon Divergence or Euclidean distance. Defaults to ``"ED"``. Use ``"JSD"`` for Jensen-Shannon Divergence.
        fitness_function (function, optional): A function by which to measure the fitness of a potential sequence. Must take a vector representing a sequence and return a float, with lower scores indicating greater fitness. Defaults to None, in which case it uses either JSD or ED.

    Returns:
        str: The generated sequence.
    '''

    if all([char in {"A", "T", "G", "C"} for char in aa_seq]):
        warn("This appears to be a DNA sequence, not an amino acid sequence. Ensure that you are passing in an amino acid sequence.")

    for target, frequencies in target_params.items():
        print(sum(frequencies.values()))
        if not isclose(sum(frequencies.values()), 1):
            raise ValueError("Target frequencies for " + str(target) + " do not sum to 1.0")

    # back translate to an initial seq
    insert = ""
    for aa in aa_seq:
        try:
            insert += Bio.Data.CodonTable.unambiguous_dna_by_id[genetic_code].back_table[aa]
        except KeyError:
            if aa == "*":
                insert += Bio.Data.CodonTable.unambiguous_dna_by_id[genetic_code].back_table[None]

    # create the genetic algorithm instance
    ga = GeneticAlgorithm(dna_to_vector(insert),
                          crossover_probability=crossover_probability,
                          maximise_fitness=False,
                          population_size=population_size,
                          mutation_probability=mutation_probability)

    # get the target values of k
    k = list(target_params.keys())

    # generate the target vector from the input dict
    target = np.array([])
    for _k in sorted([x for x in k if x != "codons"]):
        target = np.concatenate((target, [x[1] for x in sorted(target_params[_k].items(), key=lambda x: x[0])]))
    if "codons" in k:
        target = np.concatenate((target, [x[1] for x in sorted(target_params["codons"].items(), key=lambda x: x[0])]))

    def vector(seq):
        output = np.array([])
        if [k for k in k if k != "codons"]:
            output = np.concatenate((output, k_mer_frequencies(seq, [x for x in k if x != "codons"], include_missing=True, vector=True)))
        if "codons" in k:
            output = np.concatenate((output, [x[1] for x in sorted(codon_frequencies(seq).items(), key=lambda x: x[0])]))
        return output

    def fitness(individual, data):
        individual = vector_to_dna(individual)
        if mode == "JSD":
            return jensen_shannon_divergence([dit.ScalarDistribution(target / len(k)), dit.ScalarDistribution(vector(individual) / len(k))])
        if mode == "ED":
            return np.linalg.norm(target - vector(individual))
        raise Exception("Fitness mode must be JSD or ED")

    ga.fitness_function = fitness if not fitness_function else fitness_function

    synonymous_codons = _synonymous_codons(genetic_codes[genetic_code])

    def mutate(individual):
        while True:
            # choose a random codon
            codon_idx = np.random.randint(len(individual) / 6) * 6

            # figure out which codon it is
            codon = vector_to_dna(individual[codon_idx:codon_idx + 6])

            # ensure that mutations actually change the sequence
            if len(synonymous_codons[codon]) != 1:
                break

        # choose a new one at random for the AA
        new_codon = dna_to_vector(np.random.choice([x for x in synonymous_codons[codon] if x != codon]))

        # replace it in the individual
        individual[codon_idx:codon_idx + 6] = new_codon

        return individual
    ga.mutate_function = mutate

    def crossover(parent_1, parent_2):
        parent_1, parent_2 = list(parent_1), list(parent_2)
        index = random.randrange(1, len(parent_1) / 6) * 6
        child_1 = parent_1[:index] + parent_2[index:]
        child_2 = parent_2[:index] + parent_1[index:]
        return child_1, child_2
    ga.crossover_function = crossover

    def create_individual(seed_data):
        individual = vector_to_dna(seed_data)
        new = ""
        for codon in [individual[i:i + 3] for i in range(0, len(individual), 3)]:
            if len(synonymous_codons[codon]) == 1:
                new += codon
                continue
            new += np.random.choice([x for x in synonymous_codons[codon] if x != codon])

        return dna_to_vector(new)
    ga.create_individual = create_individual

    # set up for GA run
    ga.create_first_generation()
    gens_since_improvement = 0
    best_indv_fitness = ga.best_individual()[0]
    counter = 1

    # run the GA
    try:
        while gens_since_improvement < max_gens_since_improvement:
            ga.create_next_generation()
            if ga.best_individual()[0] < best_indv_fitness * (1 - improvement_rel_threshold):
                best_indv_fitness = ga.best_individual()[0]
                gens_since_improvement = 0
            else:
                gens_since_improvement += 1
            if verbose:
                print("Gen: %s\tSince Improvement: %s/%s\tFitness: %s".expandtabs(15) % (counter, gens_since_improvement, max_gens_since_improvement, ga.best_individual()[0]), end="\r")
            counter += 1
    except KeyboardInterrupt:
        print("\nStopping early...")

    if verbose:
        print()

    best_seq = vector_to_dna(ga.best_individual()[1])
    assert Seq(best_seq).translate(genetic_code) == Seq(insert).translate(genetic_code)
    return best_seq
