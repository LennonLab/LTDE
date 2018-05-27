from __future__ import division
import os, re
from biom import parse_table, load_table
#from biom import load_table
import pandas as pd
import numpy as np
import  matplotlib.pyplot as plt
import ltde_tools as lt
from collections import defaultdict


def get_soil_sites():
    df_path = lt.get_path() + '/data/EMP/EMP-clean/processed_data/spore_counts.txt'
    df = pd.read_csv(df_path, sep = '\t', index_col = 0)
    sites = df.index.values
    sites_ = [str(x) for x in sites]
    sites_str = "|".join(sites_)
    directory = lt.get_path() + '/data/EMP/EMP-public/templates'
    soil_sites = []
    for filename in os.listdir(directory):
        if filename.endswith(".txt")  & ("_prep_" not in filename):
            meta_path = os.path.join(directory, filename)
            meta = pd.read_csv(meta_path, sep = '\t', index_col = 0)
            meta.index = meta.index.map(str)
            # make sure the column is all lower case
            meta['sample_type'] = meta['sample_type'].str.lower()
            meta_soil = meta.loc[meta['sample_type'] == 'soil']
            treatment_columns = [x for x in meta.columns.values if 'treatment' in x]
            if len(treatment_columns) > 0:
                continue
            soil_sites.extend(meta.index.tolist())
    return soil_sites


def get_all_genus():
    df_path = lt.get_path() + '/data/traits/persistence.phylo.txt'
    df = pd.read_csv(df_path, sep = '\t', index_col = 0)
    genera = list(set(df.Genus.values))
    #genera = ['Terriglobus']
    genera.remove('Methanosarcina')
    genera_df_list_dict = defaultdict(list)
    to_keep = get_soil_sites()
    directory = lt.get_path() + '/data/EMP/EMP-public/processed_data'
    for filename in os.listdir(directory):
        if filename.endswith("_otu_table.biom"):
            print(filename)
            biom_path = os.path.join(directory, filename)
            table = load_table(biom_path)
            tax = table.metadata_to_dataframe('observation')
            df = table.to_dataframe()
            to_keep_biom = list(set(to_keep) & set(df.columns.tolist()))
            if len(to_keep_biom) == 0:
                continue
            df_to_keep = df.loc[:, to_keep_biom]
            N = df_to_keep.sum(0)
            for genus in genera:
                tax_genus = tax[tax['taxonomy_5'].str.contains("g__" + genus)]
                genus_otus = tax_genus.index.values
                if len(genus_otus) == 0:
                    continue
                genus_otus_ = [str(x) for x in genus_otus]
                genus_otus_str = "|".join(genus_otus_)
                df_genus = df_to_keep[df_to_keep.index.str.contains(genus_otus_str)]
                N_genus = df_genus.sum(0)
                N_merge = pd.concat([N, N_genus], axis=1)
                N_merge.columns = ['N', 'N_genus']
                N_merge = N_merge[N_merge.N_genus != 0]
                genera_df_list_dict[genus].append(N_merge)

    for g in genera:
        genus_df_list = genera_df_list_dict[g]
        if len(genus_df_list) == 0:
            continue
        genus_dfs_concat = pd.concat(genus_df_list)
        OUT_file_path = lt.get_path() + '/data/EMP/EMP-clean/processed_data/' + g + '.txt'
        genus_dfs_concat.to_csv(OUT_file_path, sep = '\t', index = True)




def plot_rel_N():
    directory = lt.get_path() + '/data/EMP/EMP-clean/processed_data'
    for filename in os.listdir(directory):
        if filename.endswith(".txt"):
            df_path = os.path.join(directory, filename)
            genus = re.split(r'[./]+', df_path)[-2]
            df = pd.read_csv(df_path, sep = '\t', index_col = 0)
            df = df[df['N'] >= 10]
            N_genus_rel = np.log10(df['N_genus'].values / df['N'].values)
            fig = plt.figure()
            plt.hist(N_genus_rel, density=True, bins=50, alpha = 0.9)
            plt.axvline(x=np.mean(N_genus_rel), c = 'k', ls = '--', lw = 3)
            plt.title(genus + ' in soil samples', fontsize = 14)
            plt.xlabel('Relative abundance, ' + r'$log_{10}$', fontsize = 12)
            plt.ylabel('Frequency', fontsize = 11)
            fig_name = lt.get_path() + '/figs/EMP/' + genus + '_rel_N.png'
            fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
            plt.close()


#get_all_genus()
plot_rel_N()
