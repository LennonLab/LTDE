from __future__ import division
import os, re
import pandas as pd
import numpy as np

mydir = os.path.expanduser("~/GitHub/LTDE/")


strains = ['ATCC13985', 'ATCC43928', 'KBS0701', 'KBS0702', 'KBS0703', 'KBS0705', \
            'KBS0706', 'KBS0707', 'KBS0710', 'KBS0711', 'KBS0712', 'KBS0714', \
            'KBS0715', 'KBS0721', 'KBS0722', 'KBS0724', 'KBS0725', 'KBS0801', \
            'KBS0802', 'KBS0812']

def module_to_KO(strain):
    kaas_directory = mydir + 'data/MAPLE/' + strain + '_MAPLE_result/KAAS'
    data = [['KEGG_Orthology', 'Pathway_ID']]
    bad_chars = '()-+,-'
    rgx = re.compile('[%s]' % bad_chars)
    for filename in os.listdir(kaas_directory):
        if filename.endswith("_matrix.txt"):
            for line in open((os.path.join(kaas_directory, filename)), 'r'):
                line_strip_split = line.strip().split()
                if len(line_strip_split) > 2 and 'M' in line_strip_split[0]:
                    if '_' in line_strip_split[0]:
                        pathway = line_strip_split[0].split('_')[0]
                    else:
                        pathway = line_strip_split[0]
                    ko_genes = line_strip_split[2:]
                    for ko_gene in ko_genes:
                        test_set_member = [bad_char for bad_char in bad_chars if bad_char in ko_gene]
                        if len(test_set_member) > 0:
                            ko_gene_clean = rgx.sub('', ko_gene)
                            ko_gene_clean_split =  ['K' + e for e in ko_gene_clean.split('K') if e]
                            for split_gene in ko_gene_clean_split:
                                if 'M' in split_gene:
                                    continue
                                data.append([split_gene, pathway])
                        else:
                            if 'K' in ko_gene:
                                data.append([ko_gene, pathway])

    df = pd.DataFrame(data[1:],columns=data[0])
    OUT_path =  mydir + 'data/MAPLE_clean/' + strain + '_MAPLE_clean/'
    if not os.path.exists(OUT_path):
        os.makedirs(OUT_path)
    OUT_file_path = OUT_path + 'KO_to_M.txt'
    df.to_csv(OUT_file_path, sep = '\t', index = False)


def clean_kaas(strain):
    IN_kaas_path = mydir + 'data/KAAS/' + strain + '_KAAS_result_ko'
    IN_kaas = pd.read_csv(IN_kaas_path, sep = '\t',
        names = ['protein_id', 'KO', 'species', 'phylum_genus', 'num'])
    IN_kaas_subset = IN_kaas.loc[IN_kaas['KO'] != 'K_NA']
    OUT_path = mydir + 'data/KAAS_clean/' + strain + '_KAAS_clean.txt'
    OUT = open(OUT_path, 'w')
    print>> OUT, 'protein_id', 'KEGG_Orthology', 'species', 'phylum', 'genus', 'num'
    count = 0
    for index, row in IN_kaas_subset.iterrows():
        KO_split =  row['KO'].split(',')
        phylum_genus_split =  row['phylum_genus'].strip().split('-')
        for KO in KO_split:
            if len(phylum_genus_split) == 1:
                # KEGG used 'Others' to for unknown genus
                genus = 'Others'
                if 'Other 'in phylum_genus_split[0]:
                    phylum = phylum_genus_split[0].replace(' ', '-')
                else:
                    phylum = phylum_genus_split[0].strip()
                    print>> OUT, row['protein_id'], KO, row['species'], \
                            phylum_genus_split[0].strip(), 'Others', str(int(row['num']))
            else:
                genus = phylum_genus_split[1].strip()
                phylum = phylum_genus_split[0].strip()
            print>> OUT, row['protein_id'], KO, row['species'], \
                    phylum, genus, str(int(row['num']))
            count += 1
    OUT.close()



def merge_maple(strain):
    IN_maple_sign_path = mydir + 'data/MAPLE/' + strain + '_MAPLE_result/' + 'module_signature.tsv'
    IN_maple_sign = pd.read_csv(IN_maple_sign_path, sep = '\t')
    IN_maple_cmplx_path = mydir + 'data/MAPLE/' + strain + '_MAPLE_result/' + 'module_complex.tsv'
    IN_maple_cmplx = pd.read_csv(IN_maple_cmplx_path, sep = '\t')
    IN_maple_pthwy_path = mydir + 'data/MAPLE/' + strain + '_MAPLE_result/' + 'module_pathway.tsv'
    IN_maple_pthwy = pd.read_csv(IN_maple_pthwy_path, sep = '\t')
    IN_maple_fxn_path = mydir + 'data/MAPLE/' + strain + '_MAPLE_result/' + 'module_function.tsv'
    IN_maple_fxn = pd.read_csv(IN_maple_fxn_path, sep = '\t')
    df_list = [IN_maple_cmplx, IN_maple_pthwy, IN_maple_sign]
    df_merged = IN_maple_fxn.append(df_list)
    # add column with pathway ID
    df_merged['Pathway_ID'] = df_merged['ID'].apply(lambda x: x.split('_')[0])
    df_merged_no_dup = df_merged.drop_duplicates(subset='Pathway_ID', keep="last")
    df_merged_no_dup = df_merged_no_dup.reset_index(drop=True)
    # median = median MCR
    OUT_path = mydir + 'data/MAPLE_clean/' + strain + '_MAPLE_clean/maple_modules.txt'
    df_merged_no_dup.to_csv(OUT_path, sep = '\t', index = False)


def merge_maple_all_strains():
    dfs = []
    maple_path = mydir + 'data/MAPLE_clean/'
    for subdir, dirs, files in os.walk(maple_path):
        for f in files:
            if f == 'maple_modules.txt':
                df = pd.read_csv(os.path.join(subdir, f), sep = '\t')
                strain = re.split(r'(_+|/+|)', subdir)[-5]
                df['Strain'] = strain
                dfs.append(df)
    dfs_concat = pd.concat(dfs)
    dfs_concat = dfs_concat.reset_index(drop=True)
    # remove rows that are less than 50% complete
    # query(coverage) = MCR % (ITR)
    #query(coverage/max) = MCR % (WC)
    #query(coverage/mode) = Q-value

    dfs_concat_050 = dfs_concat.loc[dfs_concat['query(coverage)'] >= 0.5]
    module_by_taxon = pd.crosstab(dfs_concat_050.Pathway_ID, dfs_concat_050.Strain)
    module_by_taxon_no_redundant = module_by_taxon[(module_by_taxon.T != 1).any()]
    OUT_path = mydir + 'data/module_by_taxon.txt'
    module_by_taxon_no_redundant.to_csv(OUT_path, sep = '\t', index = True)


def run_all_strains():
    for strain in strains:
        module_to_KO(strain)
        clean_kaas(strain)
        merge_maple(strain)

#run_all_strains()
#merge_maple_all_strains()
