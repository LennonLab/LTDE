from __future__ import division
import ltde_tools as lt
import glob, re, os, subprocess, math, json
import pandas as pd
import numpy as np
from collections import Counter
from Bio import SeqIO

def make_16S_fasta():
    alignments = ['KBS0710_NR_024911', 'KBS0721_NR_114994']

    def generate_16S_consenus(alignment):
        mpileup = get_path() + '/data/align/' + alignment + '.pileup'
        fasta =  get_path() + '/data/align/ltde_seqs_clean.fasta'
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


    def run_alignments():
        for align in alignments:
            generate_16S_consenus(align)

    # remove KBS0727 because it's identical to KBS0725
    #os.system(cat ~/GitHub/LTDE/data/align/ltde_seqs.fasta | awk '{if (substr($0,1) == ">KBS0727") censor=1; else if (substr($0,1,1) == ">") censor=0; if (censor==0) print $0}' > ~/GitHub/LTDE/data/align/ltde_seqs_clean.fasta)
    #run_alignments()
    # https://www.arb-silva.de/aligner/job/632523
    # NC_005042.1:353331-354795 renamed as NC_005042.1.353331-354795
    #os.system(sed -i -e 's/NC_005042.1:353331-354795/NC_005042.1.353331-354795/g' ~/GitHub/LTDE/data/align/ltde_seqs_clean.fasta)
    # ltde_neighbors_seqs.fasta uploaded to ARB and aligned
    # alignment file = arb-silva.de_2019-04-07_id632669.fasta



def get_16S_copy_number():
    genome_path = lt.get_path() + '/data/genomes/genomes_ncbi/'
    df_out = open(lt.get_path() + '/data/count_16S.txt', 'w')
    header = ['Species', 'Number_16S']
    df_out.write('\t'.join(header) + '\n')
    for subdir, dirs, files in os.walk(genome_path):
        for file in files:
            if file.endswith('.gbff'):
                strain = subdir.split('/')[-1]
                count_16S = 0
                with open(os.path.join(subdir, file), "rU") as input_handle:
                    for record in SeqIO.parse(input_handle, "genbank"):
                            for feature in record.features:
                                if feature.type == 'rRNA':
                                    if feature.qualifiers['product'][0] == '16S ribosomal RNA':
                                        count_16S+=1
                df_out.write('\t'.join([strain, str(count_16S)]) + '\n')
    df_out.close()




def clean_iRep():
    # very low coverage for these taxa
    to_remove = ['KBS0705', 'KBS0706']
    directory = os.fsencode(lt.get_path() + '/data/iRep')
    df_out = open(lt.get_path() + '/data/iRep_clean.txt', 'w')
    header = ['Sample', 'Species', 'rep' ,'iRep']
    df_out.write('\t'.join(header) + '\n')
    iRep_corrected_dict = {}
    iRep_uncorrected_dict = {}
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith('.tsv'):
            iRep_path = os.path.join(str(directory, 'utf-8'), filename)
            strain = re.split(r'[.-]+', filename)[0]
            strain_rep = re.split(r'[.]+', filename)[0]
            if strain in to_remove:
                continue
            if 'W' in strain_rep:
                continue
            #    strain_rep = 'KBS0711-5'
            #elif 'WB' in strain_rep:
            #    strain_rep = 'KBS0711-6'
            #elif 'WC' in strain_rep:
            #    strain_rep = 'KBS0711-7'
            #else:
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
    directory = os.fsencode(lt.get_path() + '/data/genomes/nanopore_hybrid_annotated_cogs')
    cog_dict = {}
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith('_reformat.txt'):
            cog_path = os.path.join(str(directory, 'utf-8'), filename)
            df = pd.read_csv(cog_path, sep = '\t')
            df_cogs = df.loc[df['source'] == 'COG_FUNCTION']
            cogs = df.accession.values
            cogs = [cog.split('!!!')[0] for cog in cogs if 'COG' in cog]
            strain = filename.split('_')[0]
            cog_dict[strain] = {}
            for cog in cogs:
                cog_dict[strain][cog] = 1

    df_cogs = pd.DataFrame.from_dict(cog_dict)
    df_cogs = df_cogs.fillna(0)
    df_cogs = df_cogs[(df_cogs.T != 1).any()]
    df_cogs = df_cogs[(df_cogs.T != 1).any()].T
    df_out = lt.get_path() + '/data/genomes/nanopore_hybrid_annotated_cogs.txt'
    df_cogs.to_csv(df_out, sep = '\t', index = True)




def merge_maple(strain):
    maple_path = lt.get_path() + '/data/genomes/genomes_ncbi_maple/'
    IN_maple_sign_path = maple_path + strain + '_MAPLE_result/' + 'module_signature.tsv'
    IN_maple_sign = pd.read_csv(IN_maple_sign_path, sep = '\t')
    IN_maple_cmplx_path = maple_path + strain + '_MAPLE_result/' + 'module_complex.tsv'
    IN_maple_cmplx = pd.read_csv(IN_maple_cmplx_path, sep = '\t')
    IN_maple_pthwy_path = maple_path + strain + '_MAPLE_result/' + 'module_pathway.tsv'
    IN_maple_pthwy = pd.read_csv(IN_maple_pthwy_path, sep = '\t')
    IN_maple_fxn_path = maple_path + strain + '_MAPLE_result/' + 'module_function.tsv'
    IN_maple_fxn = pd.read_csv(IN_maple_fxn_path, sep = '\t')
    df_list = [IN_maple_cmplx, IN_maple_pthwy, IN_maple_sign]
    df_merged = IN_maple_fxn.append(df_list)
    # add column with pathway ID
    df_merged['Pathway_ID'] = df_merged['ID'].apply(lambda x: x.split('_')[0])
    df_merged_no_dup = df_merged.drop_duplicates(subset='Pathway_ID', keep="last")
    df_merged_no_dup = df_merged_no_dup.reset_index(drop=True)
    # median = median MCR
    OUT_path = lt.get_path() + '/data/genomes/genomes_ncbi_maple_clean/' + strain + '_maple_modules.txt'
    df_merged_no_dup.to_csv(OUT_path, sep = '\t', index = False)



def merge_maple_all_strains():
    dfs = []
    maple_path = lt.get_path() + '/data/genomes/genomes_ncbi_maple_clean/'
    for filename in os.listdir(maple_path):
        if filename.endswith("_maple_modules.txt"):
            df = pd.read_csv(maple_path + filename, sep = '\t')
            strain = filename.split('_')[0]
            df['Strain'] = strain
            dfs.append(df)

    dfs_concat = pd.concat(dfs)
    dfs_concat = dfs_concat.reset_index(drop=True)
    # remove rows that are less than 50% complete
    # query(coverage) = MCR % (ITR)
    #query(coverage/max) = MCR % (WC)
    #query(coverage/mode) = Q-value
    dfs_concat_050 = dfs_concat.loc[dfs_concat['query(coverage)'] >= 0.8]
    module_by_taxon = pd.crosstab(dfs_concat_050.Pathway_ID, dfs_concat_050.Strain)
    module_by_taxon_no_redundant = module_by_taxon[(module_by_taxon.T != 1).any()]
    OUT_path = lt.get_path() + '/data/genomes/genomes_ncbi_maple.txt'
    module_by_taxon.to_csv(OUT_path, sep = '\t', index = True)


def get_assembly_coverage():
    df_out = open(lt.get_path() + '/data/genomes/assembly_coverage.txt', 'w')
    df_out.write('\t'.join(['strain', 'mean_coverage']) + '\n')
    assembly_path = lt.get_path() + '/data/genomes/nanopore_hybrid/'
    for file in os.listdir(assembly_path):
        filename = os.fsdecode(file)
        if filename.endswith('.fasta'):
            strain = filename.split('.')[0]
            print(strain)
            fa = lt.classFASTA(assembly_path+filename).readFASTA()
            fa_headers = [x[0].split('_') for x in fa]
            fa_headers = [x for x in fa_headers if int(x[3]) > 200]
            size = sum(int(x[3]) for x in fa_headers)
            weighted_mean_cov = (sum( [int(x[3]) * float(x[5]) for x in fa_headers]) / size)
            df_out.write('\t'.join([ strain, str(round(weighted_mean_cov, 3)) ]) + '\n')

    df_out.close()


def get_breseq_samples_to_ignore():
    json_path = lt.get_path() + '/data/breseq/summary/'
    to_keep = []
    for filename in os.listdir(json_path):
        if filename.endswith(".json") == False:
            continue
        with open(json_path + filename) as f:
            data = json.load(f)
            contigs = list(data['references']['reference'].keys())
            coverages = []
            for contig in contigs:
                if data['references']['reference'][contig]['length'] < 300:
                    continue
                coverages.append(data['references']['reference'][contig]['coverage_average'] )
            mean_cov = np.mean(coverages)
            if mean_cov > 20:
                to_keep.append(filename.split('.')[0])
    return to_keep





def clean_breseq():
    output_to_keep = ['INS', 'DEL', 'SNP']
    to_keep_samples = get_breseq_samples_to_ignore()
    evidence_path = lt.get_path() + '/data/breseq/output/'
    # get list of taxa to analyze
    to_keep_taxa = list(set([ x.split('-')[0] for x in to_keep_samples ]))
    for taxon in to_keep_taxa:
        taxon_samples = [ x for x in to_keep_samples if x.startswith(taxon) ]
        if len(taxon_samples) < 3:
            to_keep_taxa.remove(taxon)
    to_keep_taxa = ['KBS0715']
    for taxon in to_keep_taxa:
        taxon_sites = []
        taxon_samples = [ x for x in to_keep_samples if x.startswith(taxon) ]
        for taxon_sample in taxon_samples:
            for i, line in enumerate(open(evidence_path + taxon_sample + '.gd', 'r')):
                line_split = line.strip().split('\t')
                if line_split[0] in output_to_keep:
                    # a lot of mutations at the first base of each contig, ignore these
                    if line_split[4] == '1':
                        continue
                    taxon_sites.append( line_split[3] + '_' + str(line_split[4]))

        count_dict = Counter(taxon_sites)
        counts = list(count_dict.values())
        count_dict_4 = dict((k, v) for k, v in count_dict.items() if v >= 3)
        sites = list(count_dict_4.keys())
        #print(sites)
        for taxon_sample in taxon_samples:
            for i, line in enumerate(open(evidence_path + taxon_sample + '.gd', 'r')):
                line_split = line.strip().split('\t')
                if line_split[0] in output_to_keep:
                    if line_split[4] == '1':
                        continue
                    if line_split[3] + '_' + str(line_split[4]) in sites:
                        print(line_split)


        for i in list(range(1,5)):
            print(taxon, i, counts.count(i)/len(counts))




#clean_iRep()
#get_assembly_coverage()
#for strain in lt.strain_list():
#    merge_maple(strain)
#merge_maple_all_strains()

get_16S_copy_number()
