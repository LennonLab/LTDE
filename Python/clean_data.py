from __future__ import division
import ltde_tools as lt
import glob, re, os, subprocess, math
import pandas as pd


def make_16S_fata():

    alignments = ['KBS0710_NR_024911', 'KBS0721_NR_114994']

    def generate_16S_consenus(alignment):
        mpileup = get_path() + '/data/align/' + alignment + '.pileup'
        fasta =  get_path() + '/data/align/ltde_seqs.fasta'
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
        out_fasta.write('\n')
        out_fasta.write('>' + strain + '\n')
        split_seq = ltde_tools.split_by_n(seq_str, 60)
        for split_seq_i in split_seq:
            out_fasta.write(split_seq_i + '\n')
        out_fasta.close()



    def remove_dups():
        class_fa = ltde_tools.classFASTA(get_path()+'/data/align/ltde_neighbors/merged_neighbors.txt')
        out_fasta = open(get_path()+'/data/align/ltde_neighbors/merged_neighbors_no_dups.txt', 'w+')

        read_fa = class_fa.readFASTA()
        names = [x[0] for x in read_fa]
        def duplicates(lst, item):
            return [i for i, x in enumerate(lst) if x == item]
        dups = dict((x, duplicates(names, x)) for x in set(names) if names.count(x) > 1)
        to_remove = []
        for k, v in sorted(dups.items()):
            to_remove.extend(v[1:])
        read_fa_no_dups = [i for j, i in enumerate(read_fa) if j not in to_remove]

        for x in read_fa_no_dups:

            name = x[0]
            out_fasta.write('>' + name + '\n')
            seq = x[1]
            n = 80
            #seq_split = [seq[i:i+n] for i in range(0, len(seq), n)]
            seq_split = lt.split_by_n(seq, n)
            for seq_i in seq_split:
                out_fasta.write(seq_i + '\n')
            out_fasta.write('\n')
        out_fasta.close()


    def run_alignments():
        for align in alignments:
            generate_16S_consenus(align)

    # un-blank the commands to merge the sequence data
    #os.system('cat ~/GitHub/LTDE/data/align/ltde_neighbors/*.txt > ~/GitHub/LTDE/data/align/ltde_neighbors/merged_neighbors.txt')
    #remove_dups()
    #os.system('cat ~/GitHub/LTDE/data/align/ltde_seqs.fasta ~/GitHub/LTDE/data/align/ltde_neighbors/merged_neighbors_no_dups.txt > ~/GitHub/LTDE/data/align/ltde_neighbors_seqs.fasta')

    # https://www.arb-silva.de/aligner/job/632523
    # NC_005042.1:353331-354795 renamed as NC_005042.1.353331-354795
    os.system(sed -i -e 's/NC_005042.1:353331-354795/NC_005042.1.353331-354795/g' ~/GitHub/LTDE/data/align/ltde_neighbors_seqs.fasta)

    # ltde_neighbors_seqs.fasta uploaded to ARB and aligned
    # alignment file = arb-silva.de_2019-04-06_id632523.fasta



def get_iRep():
    directory = os.fsencode(lt.get_path() + '/data/bwa_sam')
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if (filename.endswith('C1.sam') == True) or (filename.endswith('C2.sam') == True):
            continue
        if filename.endswith('.sam'):
            strain = re.split(r'[.-]+', filename)[0]
            print(filename)
            fna_path = glob.glob(lt.get_path() + '/data/genomes/*/' + strain + '/G-Chr1.fna')[0]
            out_file = lt.get_path() + '/data/iRep/' + filename.split('.')[0]
            sam = os.path.join(str(directory, 'utf-8'), filename)
            subprocess.call(['iRep', '-f', fna_path, '-s', sam, '-o', str(out_file)])


def clean_iRep():
    to_remove = ['KBS0705', 'KBS0706']
    directory = os.fsencode(lt.get_path() + '/data/iRep')
    df_out = open(lt.get_path() + '/data/iRep_clean.txt', 'w')
    header = ['Sample', 'strain', 'rep' ,'iRep']
    df_out.write('\t'.join(header) + '\n')
    iRep_corrected_dict = {}
    iRep_uncorrected_dict = {}
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith('.tsv'):
            iRep_path = sam = os.path.join(str(directory, 'utf-8'), filename)
            strain = re.split(r'[.-]+', filename)[0]
            strain_rep = re.split(r'[.]+', filename)[0]
            if strain in to_remove:
                continue
            if 'WA' in strain_rep:
                strain_rep = 'KBS0711-5'
            elif 'WB' in strain_rep:
                strain_rep = 'KBS0711-6'
            elif 'WC' in strain_rep:
                strain_rep = 'KBS0711-7'
            else:
                strain_rep = strain_rep[:-1] + str(lt.rename_rep()[strain_rep[-1]])
            for i, line in enumerate(open(iRep_path, 'r')):
                if i == 2:
                    last_item = line.strip().split()[-1]
                    if last_item == 'n/a':
                        iRep_corrected = float('nan')
                    else:
                        iRep_corrected = float(last_item)
                    iRep_corrected_dict[strain_rep] = [iRep_corrected]
                elif i == 6:
                    iRep_uncorrected = float(line.strip().split()[-1])
                    iRep_uncorrected_dict[strain_rep] = [iRep_uncorrected]
    for key, value in iRep_corrected_dict.items():
        value.extend(iRep_uncorrected_dict[key])
    for key, value in iRep_corrected_dict.items():
        if value[1] > 11:
            continue

        if math.isnan(value[0]) == True:
            iRep = value[1]
        else:
            iRep = value[0]
        out_line = [key, key.split('-')[0], key.split('-')[1], str(iRep)]
        df_out.write('\t'.join(out_line) + '\n')


    df_out.close()


def clean_COGs():
    directory = os.fsencode(lt.get_path() + '/data/COGs')
    cog_dict = {}
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith('_hmms.hits.txt'):
            cog_path = sam = os.path.join(str(directory, 'utf-8'), filename)
            df = pd.read_csv(cog_path, sep = ',')
            strain = filename.split('_')[0]
            cog_dict[strain] = {}
            cog_list = [x for x in df.gene_hmm_id.tolist() if x != '-']
            for cog in cog_list:
                cog_dict[strain][cog] = 1

    df_cogs = pd.DataFrame.from_dict(cog_dict)
    df_cogs = df_cogs.fillna(0)
    #df_cogs = df_cogs[(df_cogs.T != 1).any()]
    #df_cogs = df_cogs[(df_cogs != 1).any()]
    #print(df_cogs)
    df_cogs = df_cogs[(df_cogs.T != 1).any()].T
    df_out = lt.get_path() + '/data/COGs/cog_by_genome.txt'
    df_cogs.to_csv(df_out, sep = '\t', index = True)



#clean_iRep()
#clean_COGs()
