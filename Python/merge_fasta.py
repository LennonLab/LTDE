from __future__ import division
import os
import skbio
from skbio import DNA

mydir = os.path.expanduser("~/GitHub/LTDE/")

fa_path = mydir + 'data/align/reseq_sanger.fasta'

fa = skbio.io.read(fa_path, format='fasta')

KBS0710_8F = DNA.read(fa_path, seq_num=1)
KBS0710_1492R = DNA.read(fa_path, seq_num=2).complement(reverse=True)
KBS0721_8F = DNA.read(fa_path, seq_num=3)
KBS0721_1492R = DNA.read(fa_path, seq_num=4).complement(reverse=True)

print(KBS0721_8F)
print(KBS0721_1492R)

#print(len(KBS0721_8F))
#print(len(KBS0721_1492R))
