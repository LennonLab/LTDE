from __future__ import division
import os

mydir = os.path.expanduser("~/GitHub/LTDE/")

alignments = ['KBS0710_NR_024911', 'KBS0721_NR_114994']


def split_by_n( seq, n ):
    """A generator to divide a sequence into chunks of n units."""
    while seq:
        yield seq[:n]
        seq = seq[n:]

def generate_16S_consenus(alignment):
    mpileup = mydir + 'data/align/' + alignment + '.pileup'
    fasta =  mydir + 'data/align/ltde_seqs.fasta'
    out_fasta = open(fasta,'a+')
    strain = alignment.split('_')[0]
    seq = []
    with open(mpileup, "r") as file:
        array = []
        for line in file:
            line_split = line.split()
            if 'KBS0710' in mpileup:
                exclude = range(436, 456) + range(1,30) + range(1519, 1523)
            elif 'KBS0721' in mpileup:
                exclude = range(162, 186) + [1, 2, 1501]

            if int(line_split[1]) in exclude:
                continue
            else:
                ref = line_split[2]
                # remove lower case
                align_site = [x.upper() for x in line_split[4]]
                coverage = len(align_site)
                A_count = align_site.count('A') / coverage
                C_count = align_site.count('C') / coverage
                G_count = align_site.count('G') / coverage
                T_count = align_site.count('T') / coverage
                nucs = ['A', 'C', 'G', 'T']
                nuc_freqs = [A_count, C_count, G_count, T_count]
                nuc_freq_max = max(nuc_freqs)
                nuc_freq_max_index = nuc_freqs.index(nuc_freq_max)
                nuc_max = nucs[nuc_freq_max_index]
                if nuc_freq_max < 0.5:
                    seq.extend(ref)
                else:
                    seq.extend(nuc_max)

    seq_str = ''.join(seq)
    print>> out_fasta, '\n'
    print>> out_fasta, '>' + strain
    split_seq = split_by_n(seq_str, 60)
    for split_seq_i in split_seq:
        print>> out_fasta, split_seq_i
    out_fasta.close()

for align in alignments:
    generate_16S_consenus(align)
