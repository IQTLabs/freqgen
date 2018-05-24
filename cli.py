import click
from Bio import SeqIO
import Bio.Data.CodonTable
import yaml

from freqgen import *

@click.group()
def cli():
    pass

@cli.command(help="Featurize a FASTA file")
@click.argument('filepath', click.Path(exists=True, dir_okay=False))
@click.option('-k', multiple=True, type=int, help="Values of k to featurize the seqs for. May be repeated.")
@click.option("-c", "--codon-usage", is_flag=True, help="Whether to include a codon frequency featurization.")
@click.option("-t", "--trans-table", type=int, default=11, help="The translation table to use. Defaults to 11, the standard genetic code.")
def featurize(filepath, k, codon_usage, trans_table):
    # get the sequences as strs
    seqs = []
    with open(filepath, "rU") as handle:
        for seq in SeqIO.parse(handle, "fasta"):
            seqs.append(str(seq.seq))

    # get the k-mer k_mer_frequencies
    for _k in k:
        print(yaml.dump({_k: k_mer_frequencies(seqs, _k)}, default_flow_style=False), end="")

    # get the codon usage frequencies
    if codon_usage:
        print(yaml.dump(dict(codons=codon_frequencies("".join(seqs), trans_table)), default_flow_style=False), end="")

@cli.command(help="Generate an amino acid sequence from FASTA")
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
    with open(filepath, "rU") as handle:
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
