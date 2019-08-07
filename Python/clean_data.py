from __future__ import division
import ltde_tools as lt
import glob, re, os, subprocess, math, json
import pandas as pd
import numpy as np
from collections import Counter
from Bio import SeqIO
import _pickle as pickle
from decimal import Decimal
import statsmodels.stats.multitest as mt
from scipy.stats import t


#np.random.seed(123456789)

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




def get_diversity_stats():
    df_out = open(lt.get_path() + '/data/breseq/genetic_diversity.txt', 'w')
    df_out.write('\t'.join(['Species', 'Sample', 'Pi', 'Theta', 'Tajimas_D', 'dN_dS_total', 'dN_dS_fixed', 'mean_freq']) + '\n')
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
    # all the diversity measures
    taxa_all = []
    n_muts_all = []
    mean_freq_list_all = []
    pi_list_all = []
    theta_list_all = []
    TD_list_all = []
    dnds_total_list_all = []
    tt_all = []
    p_value_all = []
    n_reps_all = []
    n_syn_non_muts_all = []
    for taxon in to_keep_taxa:
        print(taxon)
        effective_gene_lengths, Lsyn, Lnon, substitution_specific_synonymous_fraction = lt.calculate_synonymous_nonsynonymous_target_sizes(taxon)
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
        # list of diversity statistics
        mean_freq_list = []
        pi_list = []
        theta_list = []
        TD_list = []
        dnds_total_list = []
        n_muts_list = []
        n_syn_non_muts_list = []
        for taxon_sample in taxon_samples:
            n_c = lt.get_final_pop_size(taxon_sample)
            # get SNP identifiers
            SNP_IDs = []
            fixed_SNP_IDs = []
            for i, line in enumerate(open(lt.get_path() + '/data/breseq/output/' + taxon_sample + '.gd', 'r')):
                line_split = line.strip().split('\t')
                if line_split[0] == 'SNP':
                    if line_split[3] + '_' + line_split[4] in sites_to_remove:
                        continue
                    # fixed mutations don't count towards polymorphisms
                    if float(line_split[6].split('=')[1]) == float(1):
                        fixed_SNP_IDs.append(line_split[2])
                    else:
                        SNP_IDs.append(line_split[2])

            # go back through the file again and get the coverage info from RA lines
            freq_list = []
            n_muts = 0
            for i, line in enumerate(open(lt.get_path() + '/data/breseq/output/' + taxon_sample + '.gd', 'r')):
                line_split = line.strip().split('\t')
                if (line_split[0] == 'RA') and (line_split[1] in SNP_IDs):
                    major_cov = int(line_split[15].split('=')[1].split('/')[0]) + int(line_split[15].split('=')[1].split('/')[1])
                    minor_cov = int(line_split[18].split('=')[1].split('/')[0]) + int(line_split[18].split('=')[1].split('/')[1])
                    total_cov = int(line_split[-1].split('=')[1].split('/')[0]) + int(line_split[-1].split('=')[1].split('/')[1])
                    freq = float(line_split[20].split('=')[1])
                    freq_list.append([freq, total_cov, major_cov, minor_cov])
                    n_muts += 1

            n_muts_list.append(n_muts)

            pi = lt.get_pi(freq_list, n_c, genome_size)
            theta = lt.get_theta(freq_list, n_c, genome_size)
            mean_freq = np.mean([ float(i[0]) for i in  freq_list])
            TD = lt.get_TD(freq_list=freq_list, pi=pi, theta=theta, n_c=n_c)
            non_total = 0
            syn_total = 0
            non_fixed = 0
            syn_fixed = 0
            n_syn_non_muts = 0
            for i, line in enumerate(open(lt.get_path() + '/data/breseq/annotated/' + taxon_sample + '.gd', 'r')):
                line_split = line.strip().split('\t')
                # don't count mutations that may be ancestral
                # don't count mutations in non-coding regions or psuedoregions
                if (line_split[0] != 'SNP') or ('frequency' in line_split[6]) or (line_split[3] + '_' + line_split[4] in sites_to_remove):
                    continue
                freq = float([s for s in line_split if 'frequency=' in s][0].split('=')[1])
                n_syn_non_muts += 1
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
            n_syn_non_muts_list.append(n_syn_non_muts)
            dnds_total = ((non_total+1)/(syn_total+1))/((Lnon+1)/(Lsyn+1))
            dnds_fixed = ((non_fixed+1)/(syn_fixed+1))/((Lnon+1)/(Lsyn+1))

            mean_freq_list.append(mean_freq)
            pi_list.append(pi)
            theta_list.append(theta)
            TD_list.append(TD)
            dnds_total_list.append(dnds_total)
            df_out.write('\t'.join([ taxon, taxon_sample, str(pi), str(theta), str(TD), str(dnds_total), str(dnds_fixed), str(mean_freq)]) + '\n')

        # get taxon level stats
        n_muts_all.append(np.mean(n_muts_list))
        n_syn_non_muts_all.append(np.mean(n_syn_non_muts_list))
        mean_freq_list_all.append(np.mean(mean_freq_list))
        pi_list_all.append(np.mean(pi_list))
        theta_list_all.append(np.mean(theta_list))
        TD_list_all.append(np.mean(TD_list))
        mean_dnds_total = np.mean(dnds_total_list)
        dnds_total_list_all.append(mean_dnds_total)
        taxa_all.append(taxon)
        tt = (mean_dnds_total-1)/ (np.std(dnds_total_list) / np.sqrt(float(len(dnds_total_list))))
        p_val = t.sf(np.abs(tt), len(dnds_total_list)-1) # left-tailed one-sided t-test, so use CDF
        n_reps_all.append(len(dnds_total_list))
        tt_all.append(tt)
        p_value_all.append(p_val)

    df_out.close()

    reject, pvals_corrected, alphacSidak, alphacBonf = mt.multipletests(p_value_all, alpha=0.05, method='fdr_bh')
    print(reject, pvals_corrected, alphacSidak, alphacBonf)

    # two files, one for dnds one for rest of diversity stats
    df_out_taxa = open(lt.get_path() + '/data/breseq/genetic_diversity_taxa.txt', 'w')
    df_out_taxa.write('\t'.join(['Species', 'n_muts', 'mean_freq', 'Theta', 'Pi', 'Tajimas_D']) + '\n')
    for i in range(len(taxa_all)):
        df_out_taxa.write('\t'.join([taxa_all[i], str(n_muts_all[i]), str(mean_freq_list_all[i]), str(theta_list_all[i]), str(pi_list_all[i]), str(TD_list_all[i])]) + '\n')
    df_out_taxa.close()



    df_dNdS_taxa = open(lt.get_path() + '/data/breseq/dN_dS_taxa.txt', 'w')
    df_dNdS_taxa.write('\t'.join(['Species', 'n_reps', 'n_syn_non_muts', 'dN_dS_total', 't_stat', 'p_BH']) + '\n')
    for i in range(len(taxa_all)):
        df_dNdS_taxa.write('\t'.join([taxa_all[i], str(n_reps_all[i]), str(n_syn_non_muts_all[i]), str(dnds_total_list_all[i]), str(tt_all[i]), str(pvals_corrected[i])]) + '\n')
    df_dNdS_taxa.close()



def run_parallelism_analysis(nmin_reps=3, nmin = 3, FDR = 0.05):
    output_path = lt.get_path() + '/data/breseq/output/'
    to_keep_samples = get_breseq_samples_to_keep()
    # pass nest list with frequency, coverage of major, coverage of minor, taxon
    output_to_keep = ['INS', 'DEL', 'SNP', 'SUB']
    to_keep_samples = get_breseq_samples_to_keep()
    # get list of taxa to analyze
    to_keep_taxa = list(set([ x.split('-')[0] for x in to_keep_samples ]))
    for taxon in to_keep_taxa:
        taxon_samples = [ x for x in to_keep_samples if x.startswith(taxon) ]
        if len(taxon_samples) < nmin_reps:
            to_keep_taxa.remove(taxon)
    to_keep_taxa.sort()
    #to_keep_taxa = ['KBS0712']
    p_star_dict = {}
    G_score_list = []
    for taxon in to_keep_taxa:
        print(taxon)
        effective_gene_lengths, Lsyn, Lnon, substitution_specific_synonymous_fraction = calculate_synonymous_nonsynonymous_target_sizes(taxon)
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
        # only keep mutations that are unique to a population
        count_dict_to_remove = dict((k, v) for k, v in counts_across_reps_dict.items() if v > 1)
        sites_to_remove = list(count_dict_to_remove.keys())
        # keep insertion, deletions, and nonsynonymous SNPs
        # get size_dict
        gene_count_dict = {}
        for i, line in enumerate(open(lt.get_path() + '/data/breseq/annotated/' + taxon_sample + '.gd', 'r')):
            line_split = line.strip().split('\t')
            if (line_split[0] not in output_to_keep): #or ('frequency' in line_split[6]) or (line_split[3] + '_' + line_split[4] in sites_to_remove):
                continue
            if line_split[0] == 'SNP':
                if [s for s in line_split if 'snp_type=' in s][0].split('=')[1] == 'nonsynonymous':
                    locus_tag = [s for s in line_split if 'locus_tag=' in s][0].split('=')[1]
                    if ';' in locus_tag:
                        for locus_tag_j in locus_tag.split(';'):
                            if locus_tag_j not in gene_count_dict:
                                gene_count_dict[locus_tag_j] = 1
                            else:
                                gene_count_dict[locus_tag_j] += 1
                    else:
                        if locus_tag not in gene_count_dict:
                            gene_count_dict[locus_tag] = 1
                        else:
                            gene_count_dict[locus_tag] += 1
            else:
                if len([s for s in line_split if 'gene_position=coding' in s]) >= 1:
                    locus_tag = [s for s in line_split if 'locus_tag=' in s][0].split('=')[1]
                    if ';' in locus_tag:
                        for locus_tag_j in locus_tag.split(';'):
                            if locus_tag_j not in gene_count_dict:
                                gene_count_dict[locus_tag_j] = 1
                            else:
                                gene_count_dict[locus_tag_j] += 1

                    else:
                        if locus_tag not in gene_count_dict:
                            gene_count_dict[locus_tag] = 1
                        else:
                            gene_count_dict[locus_tag] += 1

        gene_parallelism_statistics = {}
        for gene_i, length_i in effective_gene_lengths.items():
            gene_parallelism_statistics[gene_i] = {}
            gene_parallelism_statistics[gene_i]['length'] = length_i
            gene_parallelism_statistics[gene_i]['observed'] = 0
            gene_parallelism_statistics[gene_i]['multiplicity'] = 0

        # save number of mutations and multiplicity
        for locus_tag_i, n_i in gene_count_dict.items():
            gene_parallelism_statistics[locus_tag_i]['observed'] = n_i
        # remove all genes that don't acquire mutations, excess number of zeros
        #gene_parallelism_statistics = {key:val for key, val in gene_parallelism_statistics.items() if val['observed'] > 0}
        L_mean = np.mean(list(effective_gene_lengths.values()))
        L_tot = sum(list(effective_gene_lengths.values()))
        n_tot = sum(list(gene_count_dict.values()))
        # don't include taxa with less than 20 mutations
        print("N_total = " + str(n_tot))
        if n_tot < 50:
            continue
        N_genes = len(list(effective_gene_lengths.values()))
        # go back over and calculate multiplicity
        for locus_tag_i in gene_parallelism_statistics.keys():
            gene_parallelism_statistics[locus_tag_i]['multiplicity'] = gene_parallelism_statistics[locus_tag_i]['observed'] * L_mean / effective_gene_lengths[locus_tag_i]
            gene_parallelism_statistics[locus_tag_i]['expected'] = n_tot*gene_parallelism_statistics[locus_tag_i]['length']/L_tot

        pooled_multiplicities = np.array([gene_parallelism_statistics[gene_name]['multiplicity'] for gene_name in gene_parallelism_statistics.keys() if gene_parallelism_statistics[gene_name]['multiplicity'] >=1])
        pooled_multiplicities.sort()

        pooled_tupe_multiplicities = np.array([(gene_parallelism_statistics[gene_name]['multiplicity'], gene_parallelism_statistics[gene_name]['observed']) for gene_name in gene_parallelism_statistics.keys() if gene_parallelism_statistics[gene_name]['multiplicity'] >=1])
        pooled_tupe_multiplicities = sorted(pooled_tupe_multiplicities, key=lambda x: x[0])
        pooled_tupe_multiplicities_x = [i[0] for i in pooled_tupe_multiplicities]
        pooled_tupe_multiplicities_y = [i[1] for i in pooled_tupe_multiplicities]
        pooled_tupe_multiplicities_y = [sum(pooled_tupe_multiplicities_y[i:]) / sum(pooled_tupe_multiplicities_y) for i in range(len(pooled_tupe_multiplicities_y))]

        null_multiplicity_survival = lt.NullGeneMultiplicitySurvivalFunction.from_parallelism_statistics( gene_parallelism_statistics )
        #observed_ms, observed_multiplicity_survival = lt.calculate_unnormalized_survival_from_vector(pooled_multiplicities)
        null_multiplicity_survival_copy = null_multiplicity_survival(pooled_multiplicities)
        null_multiplicity_survival_copy = [sum(null_multiplicity_survival_copy[i:]) / sum(null_multiplicity_survival_copy) for i in range(len(null_multiplicity_survival_copy)) ]
        #threshold_idx = numpy.nonzero((null_multiplicity_survival(observed_ms)*1.0/observed_multiplicity_survival)<FDR)[0][0]

        mult_survival_dict = {'Mult': pooled_multiplicities, 'Obs_fract': pooled_tupe_multiplicities_y, 'Null_fract': null_multiplicity_survival_copy}
        mult_survival_df = pd.DataFrame(mult_survival_dict)
        mult_survival_df_out = lt.get_path() + '/data/breseq/mult_survival_curves/' + taxon + '.txt'
        mult_survival_df.to_csv(mult_survival_df_out, sep = '\t', index = True)

        # get likelihood score and null test
        observed_G, pvalue = lt.calculate_total_parallelism(gene_parallelism_statistics)
        G_score_list.append((taxon, observed_G, pvalue))

        print(observed_G, pvalue)
        if pvalue >= 0.05:
            continue
        # Give each gene a p-value, get distribution
        gene_logpvalues = lt.calculate_parallelism_logpvalues(gene_parallelism_statistics)
        # remove zeros ...
        #gene_logpvalues = {key:val for key, val in gene_logpvalues.items() if val > 0}

        pooled_pvalues = []
        for gene_name in gene_logpvalues.keys():
            if gene_parallelism_statistics[gene_name]['observed']>= nmin:
                pooled_pvalues.append( gene_logpvalues[gene_name] )
        pooled_pvalues = np.array(pooled_pvalues)
        pooled_pvalues.sort()

        null_pvalue_survival = lt.NullGeneLogpSurvivalFunction.from_parallelism_statistics( gene_parallelism_statistics, nmin=nmin)
        observed_ps, observed_pvalue_survival = lt.calculate_unnormalized_survival_from_vector(pooled_pvalues, min_x=-4)
        # Pvalue version
        threshold_idx = np.nonzero((null_pvalue_survival(observed_ps)*1.0/observed_pvalue_survival)<FDR)[0][0]
        #print(null_pvalue_survival(observed_ps)*1.0/observed_pvalue_survival)
        pstar = observed_ps[threshold_idx] # lowest value where this is true
        num_significant = observed_pvalue_survival[threshold_idx]
        # make it log base 10
        logpvalues_dict = {'P_value': observed_ps/math.log(10), 'Obs_num': observed_pvalue_survival, 'Null_num': null_pvalue_survival(observed_ps)}
        logpvalues_df = pd.DataFrame(logpvalues_dict)
        logpvalues_df_out = lt.get_path() + '/data/breseq/logpvalues/' + taxon + '.txt'
        logpvalues_df.to_csv(logpvalues_df_out, sep = '\t', index = True)

        p_star_dict[taxon] = (num_significant, pstar/math.log(10))

        output_mult_gene_filename = lt.get_path() + '/data/breseq/mult_genes/' + taxon + '.txt'
        output_mult_gene = open(output_mult_gene_filename,"w")
        output_mult_gene.write(", ".join(["Gene", "Length", "Observed", "Expected", "Multiplicity", "-log10(P)"]))
        for gene_name in sorted(gene_parallelism_statistics, key=lambda x: gene_parallelism_statistics.get(x)['observed'],reverse=True):
            #if gene_name in gene_logpvalues:
            if gene_logpvalues[gene_name] >= pstar and gene_parallelism_statistics[gene_name]['observed']>=nmin:
                output_mult_gene.write("\n")
                output_mult_gene.write("%s, %0.1f, %d, %0.2f, %0.2f, %g" % (gene_name, gene_parallelism_statistics[gene_name]['length'],  gene_parallelism_statistics[gene_name]['observed'], gene_parallelism_statistics[gene_name]['expected'], gene_parallelism_statistics[gene_name]['multiplicity'], abs(gene_logpvalues[gene_name])))
        output_mult_gene.close()

    G_score_list_p_vales = [i[2] for i in G_score_list]
    reject, pvals_corrected, alphacSidak, alphacBonf = mt.multipletests(G_score_list_p_vales, alpha=0.05, method='fdr_bh')
    total_parallelism_path = lt.get_path() + '/data/breseq/total_parallelism.txt'
    total_parallelism = open(total_parallelism_path,"w")
    total_parallelism.write("\t".join(["Taxon", "G_score", "p_value", "p_value_BH"]))
    for i in range(len(pvals_corrected)):
        taxon_i = G_score_list[i][0]
        G_score_i = G_score_list[i][1]
        p_value_i = G_score_list[i][2]
        pvals_corrected_i = pvals_corrected[i]

        total_parallelism.write("\n")
        total_parallelism.write("\t".join([taxon_i, str(G_score_i), str(p_value_i), str(pvals_corrected_i)]))

    total_parallelism.close()
    with open(lt.get_path() + '/data/breseq/p_star.txt', 'wb') as file:
        file.write(pickle.dumps(p_star_dict)) # use `pickle.loads` to do the reverse


#run_parallelism_analysis()
get_diversity_stats()
