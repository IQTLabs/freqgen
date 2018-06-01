import click
from click_default_group import DefaultGroup
from Bio import SeqIO
import Bio.Data.CodonTable
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import yaml

from freqgen import k_mer_frequencies, codon_frequencies
from freqgen import generate as _generate

@click.group(cls=DefaultGroup, default='generate', default_if_no_args=True)
def freqgen():
    pass

@freqgen.command(help="Featurize a FASTA file")
@click.argument('filepath', click.Path(exists=True, dir_okay=False))
@click.option('-k', multiple=True, type=int, help="Values of k to featurize the seqs for. May be repeated.")
@click.option("-c", "--codon-usage", is_flag=True, help="Whether to include a codon frequency featurization.")
@click.option("-t", "--trans-table", type=int, default=11, help="The translation table to use. Defaults to 11, the standard genetic code.")
def featurize(filepath, k, codon_usage, trans_table):
    # get the sequences as strs
    seqs = []
    with open(filepath, "r") as handle:
        for seq in SeqIO.parse(handle, "fasta"):
            seqs.append(str(seq.seq))

    # get the k-mer k_mer_frequencies
    for _k in k:
        print(yaml.dump({_k: k_mer_frequencies(seqs, _k, include_missing=True)}, default_flow_style=False), end="")

    # get the codon usage frequencies
    if codon_usage:
        print(yaml.dump(dict(codons=codon_frequencies("".join(seqs), trans_table)), default_flow_style=False), end="")

@freqgen.command(help="Generate an amino acid sequence from FASTA")
@click.argument('filepath', click.Path(exists=True, dir_okay=False))
@click.option("--mode", type=click.Choice(["freq", "seq"]), help="Whether to use the exact AA seq or its frequencies. Defaults to freq.", default="freq")
@click.option("-t", "--trans-table", type=int, default=11, help="The translation table to use. Defaults to 11, the standard genetic code.")
@click.option("-l", "--length", type=int, help="The length of the AA sequence to generate if --mode=freq.")
@click.option("-s", "--stop-codon", is_flag=True, default=True, help="Whether to include a stop codon. Defaults to true.")
def aa(filepath, mode, trans_table, length, stop_codon):
    # translate the DNA seq, if using exact AA seq
    if mode == "seq":
        print(SeqIO.read(filepath, "fasta").seq.translate(table=trans_table))
        return

    # if using frequency mode
    seqs = []
    with open(filepath, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            try:
                aa_seq = str(record.seq.translate(table=trans_table))
            except Bio.Data.CodonTable.TranslationError:
                aa_seq = str(record.seq)
            if "*" in aa_seq:
                aa_seq.replace("*", "")
            seqs.append(aa_seq)

    aa_seq = amino_acid_seq(length, k_mer_frequencies("".join(seqs), 1))
    if stop_codon:
        aa_seq += "*"
    print(aa_seq)

@freqgen.command(help="Generate a new DNA sequence with matching features")
@click.option("-a", '--aa-seq', type=click.Path(exists=True, dir_okay=False))
@click.option("-f", '--freqs', type=click.Path(exists=True, dir_okay=False))
@click.option("-v", "--verbose", is_flag=True, default=False, help="Whether to show optimization progress. Defaults to false.")
@click.option("-i", type=int, default=50, help="How many generations to stop after no improvement. Defaults to 50.")
@click.option("-p", type=int, default=100, help="Population size. Defaults to 100.")
@click.option("-m", type=float, default=0.3, help="Mutation rate. Defaults to 0.3.")
@click.option("-t", "--trans-table", type=int, default=11, help="The translation table to use. Defaults to 11, the standard genetic code.")
@click.option("-o", '--output', type=click.Path(exists=False, dir_okay=False))
def generate(aa_seq, freqs, verbose, i, p, m, trans_table, output):
    optimized = _generate(yaml.load(open(freqs)),
                          SeqIO.read(aa_seq, "fasta").seq,
                          verbose=verbose,
                          max_gens_since_improvement=i,
                          population_size=p,
                          mutation_probability=m,
                          genetic_code=trans_table)
    if verbose:
        print("Optimized sequence:", optimized)
    if output:
        with open(output, "w+") as output_handle:
            SeqIO.write(SeqRecord(Seq(optimized), id="Optimized by Freqgen", description=""), output_handle, "fasta")
