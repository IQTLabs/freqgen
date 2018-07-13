from Bio import SeqIO
import freqgen
import json
from datetime import datetime
from statistics import mean

inserts = ["BoNT", "gfp", "insulin", "luciferase", "oxytocin", "ricin"]
targets = ["banth1", "ecoli", "ngon", "paer", "saur2", "styp2"]
k_values = [1, 2, 3, 4, 5, 6, 7]

for insert in inserts:
    for target in targets:
        for k in k_values:
            print("insert = {insert}, target = {target}, k = {k}".format(insert=insert, target=target, k=k))

            # first, get the sequence from the file
            insert_seq = SeqIO.read("example_seqs/" + insert + ".fasta", "fasta").seq

            # next, we'll get the sequences of the target genes
            target_seqs = []
            with open("example_seqs/" + target + ".heg.fasta", "r") as handle:
                for seq in SeqIO.parse(handle, "fasta"):
                    seq = str(seq.seq)
                    target_seqs.append(seq)

            target_freqs = freqgen.k_mer_frequencies(target_seqs, k)
            print("Using JSD")
            jsd_start = datetime.now()
            result_jsd = freqgen.generate({k: target_freqs}, insert_seq.translate(), verbose=True, mode="JSD")
            jsd_duration = datetime.now() - jsd_start
            result_freqs_jsd = freqgen.k_mer_frequencies(result_jsd, k)
            rel_error_jsd = []
            for key in target_freqs.keys():
                try:
                    rel_error_jsd.append(abs(target_freqs[key] - result_freqs_jsd[key]) / target_freqs[key])
                except ZeroDivisionError:
                    pass
            rel_error_jsd = mean(rel_error_jsd)

            print("Using ED")
            ed_start = datetime.now()
            result_ed = freqgen.generate({k: target_freqs}, insert_seq.translate(), verbose=True, mode="ED")
            ed_duration = datetime.now() - ed_start
            result_freqs_ed = freqgen.k_mer_frequencies(result_ed, k)
            rel_error_ed = []
            for key in target_freqs.keys():
                try:
                    rel_error_ed.append(abs(target_freqs[key] - result_freqs_ed[key]) / target_freqs[key])
                except ZeroDivisionError:
                    pass
            rel_error_ed = mean(rel_error_ed)

            output = dict(insert_freqs=freqgen.k_mer_frequencies(insert_seq, k),
                          target_freqs=target_freqs,
                          result_jsd=result_jsd,
                          result_freqs_jsd=result_freqs_jsd,
                          jsd_duration_secs=jsd_duration.total_seconds(),
                          rel_error_jsd=rel_error_jsd,
                          result_ed=result_ed,
                          result_freqs_ed=result_freqs_ed,
                          ed_duration_secs=ed_duration.total_seconds(),
                          rel_error_ed=rel_error_ed)

            with open("{insert}_to_{target}_{k}mers.json".format(insert=insert, target=target, k=k), "w+") as f:
                json.dump(output, f, indent=4)
