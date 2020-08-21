############
#
# Calculates convergence matrices for different mutation identity classes
#
# These are of the form:
#
############

import numpy
import sys, os
import ltde_tools as lt
from Bio import SeqIO

#import timecourse_utils
#import phylo_tools as pt

###################
#
# Load gene information
#
###################

#taxa = ['ATCC13985', 'KBS0702', 'KBS0707', 'KBS0715', 'KBS0721', 'KBS0722', 'KBS0724', 'KBS0801']
taxa = ['ATCC13985', 'KBS0702', 'KBS0707', 'KBS0711', 'KBS0715', 'KBS0721', 'KBS0722', 'KBS0724', 'KBS0801']
maple_types = ['signature', 'complex', 'pathway', 'function']

MCR = 0.8

kegg_dict_count = {}
maple_dict_count = {}
maple_annotation_dict = {}

treatment_count_dict = {}

for taxon in taxa:
    # make refseq to protein ID dict
    refseq_to_protein_dit = {}
    for subdir, dirs, files in os.walk(lt.get_path() + '/data/genomes/genomes_ncbi/' + taxon):
        for file in files:
            if file.endswith('.gbff'):
                with open(os.path.join(subdir, file), "rU") as input_handle:
                    for record in SeqIO.parse(input_handle, "genbank"):
                        for feature in record.features:
                            if feature.type != 'CDS':
                                continue
                            if 'incomplete' in feature.qualifiers['note'][0]:
                                continue
                            if 'frameshifted' in feature.qualifiers['note'][0]:
                                continue
                            if 'internal stop' in feature.qualifiers['note'][0]:
                                continue
                            gene_name = feature.qualifiers['locus_tag'][0]
                            if 'protein_id' in feature.qualifiers:
                                treatment_count_dict[gene_name] = feature.qualifiers['protein_id'][0]

    protein_id_kegg_dict = {}

    protein_id_kegg = open(lt.get_path() + '/data/genomes/genomes_ncbi_maple/%s_MAPLE_result/query.fst.ko' % taxon, 'r')
    # make protein ID => KEGG map
    for i, line in enumerate(protein_id_kegg):
        line = line.strip()
        items = line.split("\t")
        protein_id = items[0]
        if items[1] != 'K_NA':
            protein_id_kegg_dict[items[0]] = items[1]

    significant_genes_path=lt.get_path() + '/data/breseq/mult_genes_nonsyn_sig/%s.txt' % taxon

    if os.path.exists(significant_genes_path) == False:
        continue
    significant_genes = open(significant_genes_path, 'r')

    first_line = significant_genes.readline()
    first_line = first_line.strip()
    first_line_items = first_line.split(",")

    kegg_list = []
    # get KEGG annotation of genes with significant multiplicity
    for i, line in enumerate(significant_genes):
        line = line.strip()
        items = line.split(",")
        gene_name = items[0].strip()
        if gene_name not in treatment_count_dict:
            continue
        protein_id = treatment_count_dict[gene_name]
        if protein_id in protein_id_kegg_dict:
            kegg_list.append(protein_id_kegg_dict[protein_id])

            if protein_id_kegg_dict[protein_id] not in kegg_dict_count:
                kegg_dict_count[ protein_id_kegg_dict[protein_id]] = []

            if taxon not in kegg_dict_count[protein_id_kegg_dict[protein_id]]:
                kegg_dict_count[protein_id_kegg_dict[protein_id]].append(taxon)
    # map KEGG genes onto pathays
    # get pathways
    # go through signatures, complexes, pathways, and functions
    maple_to_keep = {}
    for maple_type in maple_types:
        maple_file = open(lt.get_path() + '/data/genomes/genomes_ncbi_maple/%s_MAPLE_result/module_%s.tsv' % (taxon, maple_type) , 'r')

        first_line_maple = maple_file.readline()
        first_line_maple = first_line_maple.strip()
        first_line_items_maple = first_line_maple.split("\t")

        for i, line in enumerate(maple_file):
            line = line.strip()
            items = line.split("\t")

            if float(items[10]) < MCR:
                continue
            # remove rows that are less than 80% complete
            # query(coverage) = MCR % (ITR)
            # query(coverage/max) = MCR % (WC)
            # query(coverage/mode) = Q-value

            maple_name = items[0].split('_')[0]
            maple_to_keep[maple_name] = {}#[items[1], items[2]]
            maple_to_keep[maple_name]['type'] = items[1]
            maple_to_keep[maple_name]['description'] = items[2]

            if maple_name not in maple_annotation_dict:
                maple_annotation_dict[maple_name] = {}
                maple_annotation_dict[maple_name]['type'] = items[1]
                maple_annotation_dict[maple_name]['description'] = items[2]

            maple_matrix = open(lt.get_path() + '/data/genomes/genomes_ncbi_maple//%s_MAPLE_result/KAAS/%s_matrix.txt' % (taxon, maple_name) , 'r')
            maple_kegg = []
            # only get kegg genes that exist in the taxon's genome
            for matrix_line_i, matrix_line in enumerate(maple_matrix):
                matrix_line = matrix_line.strip()
                matrix_items = matrix_line.split('\t')
                if 'K' not in matrix_items[0]:
                    continue
                if len(matrix_items) < 3:
                    continue
                if int(matrix_items[1]) == 0:
                    continue
                maple_kegg.append(matrix_items[0])

            maple_to_keep[maple_name]['KEGG'] = maple_kegg

    # now map kegg to pathways
    for kegg_i in kegg_list:

        for maple_i, maple_i_dict in maple_to_keep.items():

            if kegg_i in maple_i_dict['KEGG']:

                if maple_i not in maple_dict_count:
                    maple_dict_count[ maple_i] = []

                if taxon not in maple_dict_count[maple_i]:
                    maple_dict_count[maple_i].append(taxon)


convergence_table = open(lt.get_path() + '/data/breseq/parallel_pathways.txt', 'w')

# alot of ribosome complex, merge those first
convergence_table.write(", ".join(["Module ID", "Module type", "Description"] + taxa ) + '\n' )

for maple_i, maple_i_taxa in maple_dict_count.items():
    if len(maple_i_taxa) < 2:
        continue
    type = maple_annotation_dict[maple_i]['type']
    description = maple_annotation_dict[maple_i]['description']
    description = description.replace(',', ';')

    description = description.split(';')[0]
    description = description.replace(' (Pentose phosphate cycle)', '')
    description = description.replace(' (PE)', '')
    description = description.replace(' deoxyribonuleotide', '')

    if 'Ribosome' in description:
        continue

    print(maple_i, maple_i_taxa)


    taxon_presence_absence = []
    for taxon in taxa:
        if taxon in maple_i_taxa:
            taxon_presence_absence.append('+')
        else:
            taxon_presence_absence.append('-')

    convergence_table.write(", ".join([maple_i, type, description] + taxon_presence_absence) + '\n' )

convergence_table.close()
