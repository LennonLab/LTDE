from __future__ import division
import os, re
import numpy as np
import pandas as pd

mydir = os.path.expanduser("~/GitHub/LTDE/")



merged_on = ['seq_id', 'position', 'gene_list', \
    'gene_name', 'gene_position', 'gene_product', \
    'locus_tag', 'gene_strand', 'transl_table', 'reference']

to_rename = ['mutation', 'codon_ref_seq', 'codon_new_seq', 'codon_number', \
    'codon_position', 'aa_position', 'aa_ref_seq', 'aa_new_seq', \
    'frequency', 'total_cov', 'number', 'file_number', \
    'prediction', 'consensus_score', 'polymorphism_score', \
    'fisher_strand_p_value', 'ks_quality_p_value', 'bias_e_value', \
    'bias_p_value', 'reject', 'snp_type', 'type', \
    'major_base', 'minor_base', 'sample']


class cleanBreseq_annotated:

    def __init__(self, path):
        self.path = path

    def clean_value(self, row, value_name):
        value_name = value_name + '='
        data_indices = [i for i, s in enumerate(row) if '=' in s]
        gene_name_index = [i for i, s in enumerate(row) if value_name in s]
        if len(gene_name_index) > 0:
            gene_name_index = gene_name_index[0]
            gene_name_index_end = data_indices.index(gene_name_index)
            if gene_name_index_end+1 >= len(data_indices):
                gene_name = row[gene_name_index:]
                gene_name = '_'.join(gene_name)
                gene_name =  re.sub(r'[^\x00-\x7F]+','-', gene_name)
                return row[:gene_name_index] +  [gene_name]

            else:
                gene_name = row[gene_name_index:data_indices[gene_name_index_end+1]]
                gene_name = '_'.join(gene_name)
                gene_name =  re.sub(r'[^\x00-\x7F]+','-', gene_name)
                return row[:gene_name_index] +  [gene_name] + row[data_indices[gene_name_index_end+1]:]
        else:
            return row

    def variant_line(self, row, columns):
        row = self.clean_value(row, 'gene_product')
        row = self.clean_value(row, 'gene_position')
        row = self.clean_value(row, 'gene_name')
        row = self.clean_value(row, 'locus_tag')
        row = self.clean_value(row, 'between')
        if row[0] == 'SNP':
            # 7 b/c added column for insertion length (which is nan for SNPs)
            columns_6_inf = columns[7:]
            row_0_6 = row[:6]
            row_6_inf = row[6:]
            row_list = []
            values_dict = dict(item.split("=") for item in row_6_inf if '=' in item)
            for column in columns_6_inf:
                if column in values_dict:
                    row_list.append(values_dict[column])
                else:
                    row_list.append(float('nan'))
            row_list = [str(x) for x in row_list]
            # add one for point mutations
            new_row = row_0_6 + ['nan'] + row_list

            return new_row
        elif row[0] == 'SUB':
            row[6], row[5] = row[5], row[6]

            columns_7_inf = columns[7:]
            row_0_7 = row[:7]
            row_7_inf = row[7:]
            row_list = []
            values_dict = dict(item.split("=") for item in row_7_inf)
            for column in columns_7_inf:
                if column in values_dict:
                    row_list.append(values_dict[column])
                else:
                    row_list.append(float('nan'))
            row_list = [str(x) for x in row_list]
            new_row = row_0_7 + row_list
            return new_row

        elif row[0] == 'INS':
            # 7 b/c added column for insertion length (which is nan for SNPs)
            columns_6_inf = columns[7:]
            row_0_6 = row[:6]
            row_6_inf = row[6:]
            row_list = []
            values_dict = dict(item.split("=") for item in row_6_inf)
            for column in columns_6_inf:
                if column in values_dict:
                    row_list.append(values_dict[column])
                else:
                    row_list.append(float('nan'))
            row_list = [str(x) for x in row_list]
            # add length of insertion
            new_row = row_0_6 + [str(len(row_0_6[-1]))] + row_list
            return new_row

        elif row[0] == 'DEL':
            # 7 b/c added column for insertion length (which is nan for SNPs)
            columns_6_inf = columns[7:]
            row_0_6 = row[:6]
            row_6_inf = row[6:]
            row_list = []
            values_dict = dict(item.split("=") for item in row_6_inf if '=' in item)
            for column in columns_6_inf:
                if column in values_dict:
                    row_list.append(values_dict[column])
                else:
                    row_list.append(float('nan'))
            row_list = [str(x) for x in row_list]
            # add length of insertion
            new_row = row_0_6 + ['nan'] + row_list
            return new_row

    def RA_line(self, row, columns):
        row = self.clean_value(row, 'gene_product')
        row = self.clean_value(row, 'gene_position')
        row = self.clean_value(row, 'gene_name')
        row = self.clean_value(row, 'locus_tag')
        row = self.clean_value(row, 'between')
        columns = columns[8:]
        row_0_8 = row[:8]
        row_8_inf = row[8:]
        row_list = []
        values_dict = dict(item.split("=") for item in row_8_inf if '=' in item)
        for column in columns:
            if column in values_dict:
                row_list.append(values_dict[column])
            else:
                row_list.append(float('nan'))
        row_list = [str(x) for x in row_list]

        new_row = row_0_8 + row_list
        return new_row

    def MC_line(self, row, columns):
        row = self.clean_value(row, 'gene_product')
        row = self.clean_value(row, 'gene_position')
        row = self.clean_value(row, 'locus_tag')
        row = self.clean_value(row, 'gene_name')
        row = self.clean_value(row, 'gene_list')
        row = self.clean_value(row, 'html_gene_name')

        columns = columns[8:]
        row_0_8 = row[:8]
        row_8_inf = row[8:]
        row_list = []
        values_dict = dict(item.split("=") for item in row_8_inf)
        for column in columns:
            if column in values_dict:
                row_list.append(values_dict[column])
            else:
                row_list.append(float('nan'))
        row_list = [str(x) for x in row_list]

        new_row = row_0_8 + row_list
        return new_row

    def JC_line(self, row, columns):
        columns = columns[10:]
        row_0_8 = row[:10]
        row_8_inf = row[10:]
        row_list = []
        values_dict = dict(item.split("=") for item in row_8_inf)
        for column in columns:
            if column in values_dict:
                row_list.append(values_dict[column])
            else:
                row_list.append(float('nan'))
        row_list = [str(x) for x in row_list]

        new_row = row_0_8 + row_list
        return new_row


    def split_annotated(self):
        path_split = self.path.split('/')
        path = '/'.join(path_split[:-3]) + '/breseq_essentials_split/' + path_split[8:9][0]
        if not os.path.exists(path):
            os.makedirs(path)

        OUT_RA = open(path + '/evidence_RA.txt', 'w')
        OUT_MC = open(path + '/evidence_MC.txt', 'w')
        OUT_JC = open(path + '/evidence_JC.txt', 'w')
        OUT_UN = open(path + '/evidence_UN.txt', 'w')
        OUT_variants = open(path + '/evidence_variants.txt', 'w')

        columns_variants = ['type','number', 'file_number', 'seq_id', 'position', \
            'mutation', 'size', 'frequency', 'gene_list', 'gene_name', 'gene_position', \
            'gene_product', 'locus_tag', 'snp_type', 'aa_new_seq', 'aa_position', \
            'aa_ref_seq', 'codon_new_seq', 'codon_number', 'codon_position', \
            'codon_ref_seq', 'gene_strand', 'transl_table', 'insert_position', 'between']

        columns_RA = ['type','number', 'misc', 'seq_id', 'position', 'change', \
            'reference', 'sample', 'total_cov', 'new_cov', 'ref_cov', 'major_cov', \
            'minor_cov', 'major_base', 'minor_base', 'prediction', 'frequency', \
            'polymorphism_frequency', 'major_frequency', 'consensus_score', \
            'polymorphism_score', 'fisher_strand_p_value', 'ks_quality_p_value', \
            'bias_e_value', 'reject', 'snp_type', 'locus_tag', 'gene_product', \
            'gene_position', 'gene_list', 'gene_name', 'bias_p_value']

        columns_MC = ['type', 'number', 'misc', 'seq_id', 'start', 'finish', \
                'zero1', 'zero2', 'left_inside_cov', 'left_outside_cov', \
                'right_inside_cov', 'right_outside_cov', 'gene_list', 'gene_name', \
                'gene_position', 'gene_product', 'locus_tag']

        columns_JC = ['type', 'number', 'misc', 'seq_id_start', 'start', 'strand'\
            'seq_id_finish', 'finish', 'unknown1', 'unknown2' 'side_1_read_count', \
            'max_right_plus', 'coverage_minus', 'prediction', 'frequency', \
            'polymorphism_frequency', 'max_min_left_plus', 'side_2_annotate_key', \
            'alignment_overlap', 'max_left_plus', 'max_min_left', 'flanking_left', \
            'max_left', 'max_min_right_plus', 'total_non_overlap_reads', \
            'coverage_plus', 'side_2_overlap', 'neg_log10_pos_hash_p_value', \
            'max_left_minus', 'side_2_redundant', 'side_1_possible_overlap_registers', \
            'side_2_coverage', 'side_1_continuation', 'max_right_minus', \
            'side_1_redundant', 'side_2_possible_overlap_registers', \
            'unique_read_sequence', 'max_min_left_minus', \
            'new_junction_read_count', 'max_pos_hash_score', 'key', \
            'side_2_continuation', 'pos_hash_score', \
            'junction_possible_overlap_registers', 'side_1_overlap', \
            'max_min_right', 'flanking_right', 'max_right', 'side_1_coverage', \
            'side_2_read_count', 'max_min_right_minus', 'new_junction_coverage', \
            'side_1_annotate_key']

        columns_UN = ['type', 'number', 'misc', 'seq_id', 'start', 'finish']

        print('\t'.join(columns_RA), file=OUT_RA)
        print('\t'.join(columns_MC), file=OUT_MC)
        print('\t'.join(columns_JC), file=OUT_JC)
        print('\t'.join(columns_UN), file=OUT_UN)
        print('\t'.join(columns_variants), file=OUT_variants)

        #print>> OUT_RA, '\t'.join(columns_RA)
        #print>> OUT_MC, '\t'.join(columns_MC)
        #print>> OUT_JC, '\t'.join(columns_JC)
        #print>> OUT_UN, '\t'.join(columns_UN)
        #print>> OUT_variants, '\t'.join(columns_variants)


        set_test = []
        with open(self.path) as f:
            for row in f:
                row = row.split()
                if len(row) < 3:
                    continue
                if row[0] == 'DEL' or row[0] == 'SNP' or \
                    row[0] == 'SUB' or row[0] == 'INS':
                    row_clean = self.variant_line(row, columns_variants)
                    #print>> OUT_variants, '\t'.join(row_clean)
                    print('\t'.join(row_clean), file=OUT_variants)
                elif row[0] == 'RA':
                    row_clean = self.RA_line(row, columns_RA)
                    #print>> OUT_RA, '\t'.join(row_clean)
                    print('\t'.join(row_clean), file=OUT_RA)
                elif row[0] == 'MC':
                    row_clean = self.MC_line(row, columns_RA)
                elif row[0] == 'JC':
                    row_clean = self.JC_line(row, columns_RA)
                elif row[0] == 'UN':
                    #print>> OUT_UN, '\t'.join(row)
                    print('\t'.join(row), file=OUT_UN)
                else:
                    continue

        OUT_RA.close()
        OUT_MC.close()
        OUT_JC.close()
        OUT_UN.close()
        OUT_variants.close()


def get_variant_annotated(variant_type, strain):
        in_path =  mydir + 'data/breseq/breseq_essentials_split/' + strain
        variants_path = in_path + '/evidence_variants.txt'
        RA_path = in_path + '/evidence_RA.txt'
        if os.path.exists(variants_path) == True:
            IN_variants = pd.read_csv(variants_path, sep = '\t', header = 'infer')
            IN_RA = pd.read_csv(RA_path, sep = '\t', header = 'infer')
            IN_variants_SNPs = IN_variants.loc[IN_variants['type'] == variant_type]
            IN_variants_SNPs = IN_variants_SNPs.drop(['size'], axis=1)
            drop_RA = ['gene_product', 'frequency', 'gene_list', 'gene_name', \
                    'gene_position', 'locus_tag', 'snp_type', 'type', 'number']
            IN_RA = IN_RA.drop(drop_RA, axis=1)
            IN_variants_SNPs.position = IN_variants_SNPs.position.astype(str)
            IN_variants_SNPs.seq_id = IN_variants_SNPs.seq_id.astype(str)

            IN_RA.position = IN_RA.position.astype(str)
            IN_RA.seq_id = IN_RA.seq_id.astype(str)
            IN_merged = pd.merge(IN_variants_SNPs, IN_RA, how='inner',
                on=['seq_id', 'position'])
            out_path = mydir + 'data/breseq/breseq_essentials_split_clean/' + strain
            if variant_type == 'INS':
                drop_dups_on = ['type', 'seq_id', 'position']
                IN_merged = IN_merged.drop_duplicates(subset=drop_dups_on, keep="first")
            if not os.path.exists(out_path):
                os.makedirs(out_path)
            #print IN_variants_SNPs
            IN_merged.to_csv(out_path + '/' + variant_type +'.txt', sep = '\t', index = False)


def merge_variant_annotated(strain, variant_type):
    count = 0
    rootdir = mydir + 'data/breseq/breseq_essentials_split_clean/'
    for folder in os.listdir(rootdir):
        if strain not in folder:
            continue
        path = rootdir  + folder + '/' + variant_type +'.txt'
        IN = pd.read_csv(path, sep = '\t', header = 'infer')
        renamed = []
        for x in to_rename:
            x_renamed = x + '_' + folder
            IN = IN.rename(columns = {x : x_renamed})
            renamed.append(x_renamed)

        if count == 0:
            merged = IN
            frequency = 'frequency_' + folder
            merged_freq = merged[frequency]
            merged.drop(labels=[frequency], axis=1,inplace = True)
            merged.insert(len(merged.columns)-1, frequency, merged_freq)
        else:
            merged_keep = renamed + merged_on
            merged = pd.merge(merged, IN[merged_keep], \
                    how='outer', on = merged_on)
        count += 1
        test = merged.columns.tolist()
        for i, column in enumerate(merged_on):
            test.remove(column)
            test.insert(i, column)
        merged = merged.reindex_axis(test, axis=1)
    sample_freqs = [x for x in merged.columns if 'frequency_' in x]
    merged_freqs = merged[sample_freqs]
    samples = merged_freqs.shape[1]
    NoDups = merged_freqs[merged_freqs.apply(lambda x:  x.isnull().sum() == samples - 1, 1)]
    NoDups_index = NoDups.index.values
    merged_unique = merged.ix[NoDups_index]
    if merged_unique.shape[0] == 0:
        return

    OUT_path = mydir + 'data/breseq/breseq_essentials_split_clean_split/' + strain
    if not os.path.exists(OUT_path):
        os.makedirs(OUT_path)

    reps = [x.split('_')[1] for x in merged_unique.columns if 'frequency_' in x]
    for rep in reps:
        to_rename_rep = [x + '_' + rep for x in to_rename]
        to_slice = merged_on + to_rename_rep
        merged_unique_rep = merged_unique[to_slice]
        set_p =  list(set(merged_unique_rep['bias_p_value_' + rep].values))
        if (len(set_p)) == 1 and (np.isnan(set_p[0]) == True):
            continue

        merged_unique_rep = merged_unique_rep[np.isfinite(merged_unique_rep['bias_p_value_' + rep])]
        OUT_path_rep = OUT_path + '/' + rep + '_' + variant_type +'.txt'
        merged_unique_rep.to_csv(OUT_path_rep, sep = '\t', index = False)


    #OUT_name = OUT_path + '/' + variant_type + '.txt'
    #merged_unique.to_csv(OUT_name, sep = '\t', index = False)


def run_everything():
    strains = []
    variants = ['SNP','INS','DEL']
    rootdir = mydir + 'data/breseq/breseq_essentials'
    for filename in os.listdir(rootdir):
        if filename == '.DS_Store':
            continue
        #annotated_path = rootdir + '/' + filename + '/annotated.gd'
        #cleanBreseq_annotated(annotated_path).split_annotated()
        #for variant in variants:
        #    get_variant_annotated(variant, strain = filename)
        strains.append(filename.split('-')[0])

    strains = list(set(strains))
    for strain in strains:
        print(strain)
        for variant in variants:
            merge_variant_annotated(strain, variant)
    # split the merged file into reps
    #if unique_mutations == True:
    #    #print "unique_mutations"
    #    get_unique_mutations(day, strain, variant_type)




run_everything()
