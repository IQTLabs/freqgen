from pyeasyga.pyeasyga import GeneticAlgorithm
from Bio import SeqIO
from freqgen import *

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
	dna = np.array([dna[i:i+2] for i in range(0, len(dna), 2)])
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
    return {codon : codons_for_amino_acid[genetic_code_dict[codon]] for codon in genetic_code_dict.keys()}

def CUB_vector(seq, genetic_code=11):
    '''returns a vector containing codon frequency'''

    frequencies = codon_frequencies(seq, genetic_code)

    # convert into vector from dict with fixed order
    frequencies = [x[1] for x in sorted(list(frequencies.items()), key=lambda x: x[0])]
    assert len(frequencies) == 64

    frequencies = np.array(frequencies)
    assert np.isclose(sum(frequencies), 1)
    return frequencies

def vector(seq):
    return np.concatenate((CUB_vector(seq),
                           # k_mer_frequencies(seq, 1, include_missing=True, vector=True),
                           # k_mer_frequencies(seq, 2, include_missing=True, vector=True),
                           # k_mer_frequencies(seq, 4, include_missing=True, vector=True),
                           # k_mer_frequencies(seq, 5, include_missing=True, vector=True),
                           k_mer_frequencies(seq, 2, include_missing=True, vector=True)
                           ))

# get target str
target = ""
for record in SeqIO.parse("ecoli.heg.fasta", "fasta"):
    target += str(record.seq)

# get insert str
insert = str(SeqIO.read("gfp.fasta", "fasta").seq)

print(f"Starting sequence: {insert}")

ga = GeneticAlgorithm(dna_to_vector(insert), crossover_probability=0, maximise_fitness=False, population_size=100, mutation_probability=0.3)

k = [3]
insert_freqs = vector(insert) # k_mer_frequencies(insert, k, include_missing=True, vector=True)
target = vector(target) # k_mer_frequencies(target, k, include_missing=True, vector=True)
print(f"Starting freqs: {insert_freqs}\nTarget freqs: {target}")
synonymous_codons = _synonymous_codons(genetic_codes[11])

def fitness(individual, data):
    individual = vector_to_dna(individual)
    fitness = np.linalg.norm(target - vector(individual)) # k_mer_frequencies(individual, k, include_missing=True, vector=True))
    return fitness
ga.fitness_function = fitness

def mutate(individual):

    while True:
        # choose a random codon
        codon_idx = np.random.randint(len(individual) / 6) * 6

        # figure out which codon it is
        codon = vector_to_dna(individual[codon_idx:codon_idx+6])

        # ensure that mutations actually change the sequence
        if len(synonymous_codons[codon]) != 1:
            break

    # choose a new one at random for the AA
    new_codon = dna_to_vector(np.random.choice([x for x in synonymous_codons[codon] if x != codon]))

    # replace it in the individual
    individual[codon_idx:codon_idx+6] = new_codon

    return individual
ga.mutate_function = mutate

def create_individual(seed_data):
    individual = vector_to_dna(seed_data)
    new = ""
    for codon in [individual[i:i+3] for i in range(0, len(individual), 3)]:
        if len(synonymous_codons[codon]) == 1:
            new += codon
            continue
        new += np.random.choice([x for x in synonymous_codons[codon] if x != codon])

    return dna_to_vector(new)
ga.create_individual = create_individual

ga.create_first_generation()

gens_since_improvement = 0
best_indv_fitness = ga.best_individual()[0]
counter = 1

while gens_since_improvement < 50:
    ga.create_next_generation()
    if ga.best_individual()[0] < best_indv_fitness:
        best_indv_fitness = ga.best_individual()[0]
        gens_since_improvement = 0
    else:
        gens_since_improvement += 1
    print(f"Gen: {counter}\tSince Improvement: {gens_since_improvement}\tFitness: {best_indv_fitness}".expandtabs(15))
    counter += 1

best_seq = vector_to_dna(ga.best_individual()[1])
print(f"Best sequence: {best_seq}")
best_freqs = vector(best_seq) # k_mer_frequencies(best_seq, k, include_missing=True, vector=True)

assert translate(best_seq) == translate(insert)

print("plotting!")
import numpy as np
import matplotlib.pyplot as plt

print("generating rectangles!")
N = len(best_freqs)
ind = np.arange(N)  # the x locations for the groups
width = 0.2       # the width of the bars
fig, ax = plt.subplots()
rects1 = ax.bar(ind, insert_freqs, width, color='b')
rects2 = ax.bar(ind + width, target, width, color='g')
rects3 = ax.bar(ind + 2*width, best_freqs, width, color='r')

print("fixing axes!")
# add some text for labels, title and axes ticks
ax.set_xticks(ind + 2*width / 2)
ax.set_xticklabels([f"{i}" for i in range(len(target))])
ax.legend((rects1[0], rects2[0], rects3[0]), ('Original', 'Target', "Optimized"))
ax.set_xlabel("k")

print("plotting!")
plt.show()
