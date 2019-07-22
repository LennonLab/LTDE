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
    df_out.write('\t'.join(['Species', 'mean_coverage']) + '\n')
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


def get_breseq_samples_to_keep():
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



def calculate_synonymous_nonsynonymous_target_sizes(taxon):
    #lines 308-334 from GitHub repo benjaminhgood/LTEE-metagenomic
    # and is under a GPL v2 license
    codon_table = { 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'CGT': 'R',
                    'CGC': 'R', 'CGA':'R', 'CGG':'R', 'AGA':'R', 'AGG':'R',
                    'AAT':'N', 'AAC':'N', 'GAT':'D', 'GAC':'D', 'TGT':'C',
                    'TGC':'D', 'CAA':'Q', 'CAG':'Q', 'GAA':'E', 'GAG':'E',
                    'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G', 'CAT':'H',
                    'CAC':'H', 'ATT':'I', 'ATC':'I', 'ATA':'I', 'TTA':'L',
                    'TTG':'L', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
                    'AAA':'K', 'AAG':'K', 'ATG':'M', 'TTT':'F', 'TTC':'F',
                    'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'TCT':'S',
                    'TCC':'S', 'TCA':'S', 'TCG':'S', 'AGT':'S', 'AGC':'S',
                    'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T', 'TGG':'W',
                    'TAT':'Y', 'TAC':'Y', 'GTT':'V', 'GTC':'V', 'GTA':'V',
                    'GTG':'V', 'TAA':'!', 'TGA':'!', 'TAG':'!' }

    # calculate number of synonymous opportunities for each codon
    codon_synonymous_opportunity_table = {}
    for codon in codon_table.keys():
        codon_synonymous_opportunity_table[codon] = {}
        for i in range(0,3):
            codon_synonymous_opportunity_table[codon][i] = -1 # since G->G is by definition synonymous, but we don't want to count it
            codon_list = list(codon)
            for base in ['A','C','T','G']:
                codon_list[i]=base
                new_codon = "".join(codon_list)
                if codon_table[codon]==codon_table[new_codon]:
                    # synonymous!
                    codon_synonymous_opportunity_table[codon][i]+=1

    codon_synonymous_substitution_table = {}
    codon_nonsynonymous_substitution_table = {}
    for codon in codon_table.keys():
        codon_synonymous_substitution_table[codon] = [[],[],[]]
        codon_nonsynonymous_substitution_table[codon] = [[],[],[]]

        for i in range(0,3):
            reference_base = codon[i]

            codon_list = list(codon)
            for derived_base in ['A','C','T','G']:
                if derived_base==reference_base:
                    continue
                substitution = '%s->%s' % (reference_base, derived_base)
                codon_list[i]=derived_base
                new_codon = "".join(codon_list)
                if codon_table[codon]==codon_table[new_codon]:
                    # synonymous!
                    codon_synonymous_substitution_table[codon][i].append(substitution)
                else:
                    codon_nonsynonymous_substitution_table[codon][i].append(substitution)
    bases = set(['A','C','T','G'])
    substitutions = []
    for b1 in bases:
        for b2 in bases:
            if b2==b1:
                continue
            substitutions.append( '%s->%s' % (b1,b2) )
    substitution_specific_synonymous_sites = {substitution: 0 for substitution in substitutions}
    substitution_specific_nonsynonymous_sites = {substitution: 0 for substitution in substitutions}
    effective_gene_synonymous_sites = {}
    effective_gene_nonsynonymous_sites = {}
    gene_length_map = {}
    genome_path = lt.get_path() + '/data/genomes/genomes_ncbi/' + taxon
    for subdir, dirs, files in os.walk(genome_path):
        for file in files:
            if file.endswith('.gbff'):
                with open(os.path.join(subdir, file), "rU") as input_handle:
                    for record in SeqIO.parse(input_handle, "genbank"):
                        for feature in record.features:
                            if feature.type != 'CDS':
                                continue
                            #print(feature)
                            if 'incomplete' in feature.qualifiers['note'][0]:
                                continue
                            if 'frameshifted' in feature.qualifiers['note'][0]:
                                continue
                            if 'internal stop' in feature.qualifiers['note'][0]:
                                continue
                            gene_name = feature.qualifiers['locus_tag'][0]

                            if gene_name not in effective_gene_synonymous_sites:
                                effective_gene_synonymous_sites[gene_name]=0
                                effective_gene_nonsynonymous_sites[gene_name]=0
                            aa_str = str(feature.qualifiers['translation'][0])
                            nuc_str = str(feature.location.extract(record).seq[:-3])
                            gene_length_map[gene_name] = len(nuc_str)
                            for position in range(len(nuc_str)):
                                codon_start = int(position/3)*3
                                codon = nuc_str[codon_start:codon_start+3]
                                if len(codon) <3:
                                    continue
                                position_in_codon = position%3

                                effective_gene_synonymous_sites[gene_name] += codon_synonymous_opportunity_table[codon][position_in_codon]/3.0
                                effective_gene_nonsynonymous_sites[gene_name] += 1-codon_synonymous_opportunity_table[codon][position_in_codon]/3.0

                                for substitution in codon_synonymous_substitution_table[codon][position_in_codon]:
                                    substitution_specific_synonymous_sites[substitution] += 1

                                for substitution in codon_nonsynonymous_substitution_table[codon][position_in_codon]:
                                    substitution_specific_nonsynonymous_sites[substitution] += 1

    substitution_specific_synonymous_fraction = {substitution: substitution_specific_synonymous_sites[substitution]*1.0/(substitution_specific_synonymous_sites[substitution]+substitution_specific_nonsynonymous_sites[substitution]) for substitution in substitution_specific_synonymous_sites.keys()}
    effective_gene_lengths = {gene_name: gene_length_map[gene_name]-effective_gene_synonymous_sites[gene_name] for gene_name in gene_length_map.keys()}
    effective_gene_lengths['synonymous'] = sum([effective_gene_synonymous_sites[gene_name] for gene_name in gene_length_map.keys()])
    effective_gene_lengths['nonsynonymous'] = sum([effective_gene_nonsynonymous_sites[gene_name] for gene_name in gene_length_map.keys()])
    effective_gene_lengths['noncoding'] = lt.get_genome_size_dict()[taxon] - effective_gene_lengths['synonymous']-effective_gene_lengths['nonsynonymous']

    return effective_gene_lengths['synonymous'], effective_gene_lengths['nonsynonymous'], substitution_specific_synonymous_fraction



def get_diversity_stats():
    df_out = open(lt.get_path() + '/data/breseq/genetic_diversity.txt', 'w')
    df_out.write('\t'.join(['Species', 'Sample', 'Pi', 'Theta', 'Tajimas_D', 'dN_dS_total', 'dN_dS_fixed']) + '\n')
    # pass nest list with frequency, coverage of major, coverage of minor, taxon
    output_to_keep = ['INS', 'DEL', 'SNP']
    to_keep_samples = get_breseq_samples_to_keep()
    output_path = lt.get_path() + '/data/breseq/output/'
    # get list of taxa to analyze
    to_keep_taxa = list(set([ x.split('-')[0] for x in to_keep_samples ]))
    for taxon in to_keep_taxa:
        taxon_samples = [ x for x in to_keep_samples if x.startswith(taxon) ]
        if len(taxon_samples) < 3:
            to_keep_taxa.remove(taxon)

    for taxon in to_keep_taxa:
        Lsyn, Lnon, substitution_specific_synonymous_fraction = calculate_synonymous_nonsynonymous_target_sizes(taxon)
        taxon_sites = []
        taxon_samples = [ x for x in to_keep_samples if x.startswith(taxon) ]
        for taxon_sample in taxon_samples:
            for i, line in enumerate(open(output_path + taxon_sample + '.gd', 'r')):
                line_split = line.strip().split('\t')
                if line_split[0] in output_to_keep:
                    # a lot of mutations at the first base of each contig, ignore these
                    if line_split[4] == '1':
                        continue
                    taxon_sites.append( line_split[3] + '_' + str(line_split[4]))

        counts_across_reps_dict = Counter(taxon_sites)
        counts_across_reps = list(counts_across_reps_dict.values())
        count_dict_to_remove = dict((k, v) for k, v in counts_across_reps_dict.items() if v > 2)
        sites_to_remove = list(count_dict_to_remove.keys())
        genome_size = lt.get_genome_size_dict()[taxon]
        taxon_sites = []
        taxon_samples = [ x for x in to_keep_samples if x.startswith(taxon) ]
        for taxon_sample in taxon_samples:
            n_c = lt.get_final_pop_size(taxon_sample)
            # get SNP identifiers
            SNP_IDs = []
            fixed_SNP_IDs = []
            for i, line in enumerate(open(lt.get_path() + '/data/breseq/output/' + taxon_sample + '.gd', 'r')):
                line_split = line.strip().split('\t')
                if line_split[0] == 'SNP':
                    #print(line_split[3] + '_' + line_split[4] in sites_to_remove)
                    if line_split[3] + '_' + line_split[4] in sites_to_remove:
                        continue
                    # fixed mutations don't count towards polymorphisms
                    if float(line_split[6].split('=')[1]) == float(1):
                        fixed_SNP_IDs.append(line_split[2])
                    else:
                        SNP_IDs.append(line_split[2])

            # go back through the file again and get the coverage info from RA lines
            freq_list = []
            for i, line in enumerate(open(lt.get_path() + '/data/breseq/output/' + taxon_sample + '.gd', 'r')):
                line_split = line.strip().split('\t')
                if (line_split[0] == 'RA') and (line_split[1] in SNP_IDs):
                    major_cov = int(line_split[15].split('=')[1].split('/')[0]) + int(line_split[15].split('=')[1].split('/')[1])
                    minor_cov = int(line_split[18].split('=')[1].split('/')[0]) + int(line_split[18].split('=')[1].split('/')[1])
                    total_cov = int(line_split[-1].split('=')[1].split('/')[0]) + int(line_split[-1].split('=')[1].split('/')[1])
                    freq = float(line_split[16].split('=')[1])
                    freq_list.append([freq, total_cov, major_cov, minor_cov])
            pi = lt.get_pi(freq_list, n_c, genome_size)
            theta = lt.get_theta(freq_list, n_c, genome_size)
            TD = lt.get_TD(freq_list=freq_list, pi=pi, theta=theta, n_c=n_c)
            non_total = 0
            syn_total = 0
            non_fixed = 0
            syn_fixed = 0
            for i, line in enumerate(open(lt.get_path() + '/data/breseq/annotated/' + taxon_sample + '.gd', 'r')):
                line_split = line.strip().split('\t')
                # don't count mutations that may be ancestral
                # don't count mutations in non-coding regions or psuedoregions
                if (line_split[0] != 'SNP') or ('frequency' in line_split[6]) or (line_split[3] + '_' + line_split[4] in sites_to_remove):
                    continue
                freq = float([s for s in line_split if 'frequency=' in s][0].split('=')[1])
                if freq == float(1):
                    if line_split[6].split('=')[1] == line_split[8].split('=')[1]:
                        syn_fixed += 1
                    else:
                        non_fixed += 1
                if line_split[6].split('=')[1] == line_split[8].split('=')[1]:
                    syn_total += 1
                else:
                    non_total += 1

            # add psuedocount of 1
            dnds_total = ((non_total+1)/(syn_total+1))/((Lnon+1)/(Lsyn+1))
            dnds_fixed = ((non_fixed+1)/(syn_fixed+1))/((Lnon+1)/(Lsyn+1))
            print(taxon_sample, dnds_total, dnds_fixed)
            df_out.write('\t'.join([ taxon, taxon_sample, str(pi), str(theta), str(TD), str(dnds_total), str(dnds_fixed)]) + '\n')

    df_out.close()



#def run_parallelism_analysis():


#calculate_synonymous_nonsynonymous_target_sizes('KBS0702')



get_diversity_stats()
