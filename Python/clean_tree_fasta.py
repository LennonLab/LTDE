from __future__ import division
import os, re
from Bio import Phylo

mydir = os.path.expanduser("~/GitHub/LTDE/")


to_replace = [('_delta/epsilon_', '_delta_epsilon_'),
                ('_Bacteroidetes/Chlorobi_', '_Bacteroidetes_Chlorobi_'),
                ('H2/3', 'H2_3'),
                ('_Chlamydiae/', '_Chlamydiae_'),
                (',_DSM_', '__DSM_'),
                ('_CP_', '_CP-'),
                ('CP_KSB1', 'CP-KSB1'),
                ('Bacteria_CPR_Parcubacteria_OD1_ALUMROCK_MS4_OD1_33_19', 'Bacteria_CP-OD1_ALUMROCK_MS4_OD1_33_19'),
                ('Chlamydiae/Verrucomicrobia', 'Chlamydiae_Verrucomicrobia'),
                ("_'Patoc_1_Paris'", '__Patoc_1_Paris_'),
                ("_baker's_yeast", '_baker_s_yeast'),
                ('_Fibrobacteres/Acidobacteria_', '_Fibrobacteres_Acidobacteria_'),
                ('_DL31_:_DL31', '_DL31___DL31'),
                ('_Chlamydia/Chlamydophila_', '_Chlamydia_Chlamydophila_'),
                ('_S26/3', '_S26_3'),
                ('_JW/YL-', '_JW_YL-'), ('_21/22', '_21_22'), ('PHE/MN1', 'PHE_MN1'),
                ('Bacteria_CP-TM6-candidate_division_TM6_bacterium_JCVI_TM6SC1__Phylum_TM6', 'Bacteria_CP_TM6-candidate_division_TM6_bacterium_JCVI_TM6SC1__Phylum_TM6'),
                ('Archaea_Euryarchaeota_Halobacteria_Haloferacales_Haloferacaceae_Halohasta_litchfieldiae_tADL_:_Hlit_replicon_OLD_NAME:halophilic_archaeon_sp._True-ADL', 'Archaea_Euryarchaeota_Halobacteria_Haloferacales_Haloferacaceae_Halohasta_litchfieldiae_tADL___Hlit_replicon_OLD_NAME_halophilic_archaeon_sp._True-ADL'),
                ('Archaea_CP-YNPFFA_candidate_division_YNPFFA_archaeon_SCGC_AAA471-B05_Combined_Assembly_YNPFFA_1__YNPFFA', 'Archaea_CP_YNPFFA_candidate_division_YNPFFA_archaeon_SCGC_AAA471-B05_Combined_Assembly_YNPFFA_1__YNPFFA'),
                ('Archaea_Euryarchaeota_Halobacteria_halophilic_archaeon_sp._DL31___DL31_replicon_1_OLD_NAME:_halophilic_archaeon_sp._DL31', 'Archaea_Euryarchaeota_Halobacteria_halophilic_archaeon_sp._DL31___DL31_replicon_1_OLD_NAME__halophilic_archaeon_sp._DL31'),
                ('Bacteria_CP-TM6-GWF2_TM6_33_332', 'Bacteria_CP_TM6-GWF2_TM6_33_332'),
                ('Bacteria_CP-TM6-gwf2_TM6_43_87', 'Bacteria_CP_TM6-gwf2_TM6_43_87'),
                ('Bacteria_Proteobacteria_Acidithiobacillia_Acidithiobacillales_Acidithiobacillaceae_Acidithiobacillus_ferrivorans_SS3', 'Bacteria_Proteobacteria_Acidithiobacillia_Acidithiobacillales_Acidithiobacillaceae_Acidithiobacillus_Acidithiobacillus_ferrivorans_SS3'),
                ('Archaea_CP-YNPFFA_candidate_division_YNPFFA_archaeon_SCGC_AAA471-B23_GBS-N_001_5', 'Archaea_CP_YNPFFA_candidate_division_YNPFFA_archaeon_SCGC_AAA471-B23_GBS-N_001_5'),
                ('Bacteria_CP-RIF31_RBG_16_RIF31_55_9', 'Bacteria_CP_RIF31_RBG_16_RIF31_55_9'),
                ('Bacteria_CP-TM6-GWE2_TM6_42_60', 'Bacteria_CP_TM6-GWE2_TM6_42_60'),
                ('Bacteria_CPR_CPR2-GWD2_CPR2_39_7', 'Bacteria_CPR2_GWD2_CPR2_39_7'),
                ('Bacteria_unclassified_Bacteria_Thermobaculum_Thermobaculum_terrenum_ATCC_BAA-798', 'Bacteria_Chloroflexi_Thermobaculum_Thermobaculum_terrenum_ATCC_BAA-798'),
                ('Bacteria_bacterium_JKG1', 'Bacteria_Chloroflexi_bacterium_JKG1'),
                ('Bacteria_CP-TM6-GWE2_TM6_36_25', 'Bacteria_CP_TM6-GWE2_TM6_36_25'),
                ('Bacteria_CP-TM7_candidate_division_TM7_genomosp_GTL1', 'Bacteria_CP_TM7_candidate_division_TM7_genomosp_GTL1'),
                ('Bacteria_Unclassified_bacteria_CG1_02_FULL_CP-38_46', 'Bacteria_Unclassified_bacteria_CG_Elusi_02'),
                ('Bacteria_CP-RIF32_BJP_IG2103_RIF32_50_23', 'Bacteria_CP-ACD39_BJP_IG2103_ACD39_50_23'),
                ('Bacteria_CPR_PeregrinibacteriaPER-GWB1_PER_54_5', 'Bacteria_CPR_PER-ii_GWB1_PER_54_5'),
                ('Bacteria_Proteobacteria_Alphaproteobacteria_Rhizobiales_Rhizobiaceae_Rhizobium/Agrobacterium_group_Rhizobium_leguminosarum_bv._trifolii_CB782', 'Bacteria_Proteobacteria_Alphaproteobacteria_Rhizobiales_Rhizobiaceae_Rhizobium_Agrobacterium_group_Rhizobium_leguminosarum_bv._trifolii_CB782'),
                ('Archaea_CP-YNPFFA_candidate_division_YNPFFA_archaeon_SCGC_AAA471-O08_GBS-N_001_29', 'Archaea_CP_YNPFFA_candidate_division_YNPFFA_archaeon_SCGC_AAA471-O08_GBS-N_001_29'),
                ('Bacteria_CP-RIF30_GWF2_RIF30_38_17', 'Bacteria_CP_RIF30_GWF2_RIF30_38_17'),
                ('Archaea_Euryarchaeota_Halobacteria_Halobacteriales_Halobacteriaceae_Natrinema_pellirubrum_157,_JCM_10476', 'Archaea_Euryarchaeota_Halobacteria_Halobacteriales_Halobacteriaceae_Natrinema_pellirubrum_157__JCM_10476'),
                ("Bacteria_Tenericutes_Mollicutes_Acholeplasmatales_Acholeplasmataceae_Candidatus_Phytoplasma_aster_yellows_witches'-broom_AY-WB", 'Bacteria_Tenericutes_Mollicutes_Acholeplasmatales_Acholeplasmataceae_Candidatus_Phytoplasma_aster_yellows_witches_-broom_AY-WB'),
                ('Bacteria_CP-TM6-GWF2_TM6_30_66', 'Bacteria_CP_TM6-GWF2_TM6_30_66'),
                ('Bacteria_Bacteroidetes_Chlorobi_group_Chlorobi_Chlorobia_Chlorobiales_Chlorobiaceae_Chlorobium/Pelodictyon_group_Pelodictyon_luteolum_DSM_273', 'Bacteria_Bacteroidetes_Chlorobi_group_Chlorobi_Chlorobia_Chlorobiales_Chlorobiaceae_Chlorobium_Pelodictyon_group_Pelodictyon_luteolum_DSM_273'),
                ('Bacteria_CP-Tectomicrobia_Candidatus_Entotheonella_sp._TSY1', 'Bacteria_CP_Tectomicrobia_Candidatus_Entotheonella_sp._TSY1'),
                ('Bacteria_CP-Tectomicrobia_Candidatus_Entotheonella_sp._TSY2', 'Bacteria_CP_Tectomicrobia_Candidatus_Entotheonella_sp._TSY2'),
                ('Archaea_CP-YNPFFA_candidate_division_YNPFFA_archaeon_SCGC_AAA471-C03_GBS-N_001_6', 'Archaea_CP_YNPFFA_candidate_division_YNPFFA_archaeon_SCGC_AAA471-C03_GBS-N_001_6'),
                ('Bacteria_CP-TM6-GWF2_TM6_32_72', 'Bacteria_CP_TM6-GWF2_TM6_32_72'),
                ('Bacteria_CP-TM6-GWE2_TM6_41_16', 'Bacteria_CP_TM6-GWE2_TM6_41_16'),
                ('Bacteria_Proteobacteria_delta_epsilon_subdivisions_Deltaproteobacteria_Myxococcales_Cystobacterineae_Myxococcaceae_Anaeromyxobacter_dehalogenans_2CP_C', 'Bacteria_Proteobacteria_delta_epsilon_subdivisions_Deltaproteobacteria_Myxococcales_Cystobacterineae_Myxococcaceae_Anaeromyxobacter_dehalogenans_2CP-C'),
                ('Bacteria_CP-TG3_candidate_division_TG3_bacterium_ACht1', 'Bacteria_CP_TG3_candidate_division_TG3_bacterium_ACht1'),
                ('Archaea_Diapherotrites_candidate_division_DUSEL3_archaeon_SCGC_AAA011-E11_DUSEL_001_206', 'Archaea_DPANN_Diapherotrites_candidate_division_DUSEL3_archaeon_SCGC_AAA011-E11_DUSEL_001_206'),
                ('Archaea_CP-DUSEL1_candidate_division_DUSEL1_archaeon_SAG-AAA011-D5_Dusel_001_270', 'Archaea_CP_DUSEL1_candidate_division_DUSEL1_archaeon_SAG-AAA011-D5_Dusel_001_270'),
                ('Bacteria_CP-Zixibacteria_BJP_IG2158_Zixibacteria_48_33', 'Bacteria_CP_Zixibacteria_BJP_IG2158_Zixibacteria_48_33'),
                ('Bacteria_CPR_Dojkabacteria_WS6_BJP_IG2158_WS6_35_49', 'Bacteria_CP-WS6_BJP_IG2158_WS6_35_49'),
                ('Bacteria_CP-TM6-GWF2_TM6_37_49', 'Bacteria_CP_TM6-GWF2_TM6_37_49')]


class classFASTA:

    def __init__(self, fileFASTA):
        self.fileFASTA = fileFASTA

    def readFASTA(self):
        '''Checks for fasta by file extension'''
        file_lower = self.fileFASTA.lower()
        '''Check for three most common fasta file extensions'''
        if file_lower.endswith('.txt') or file_lower.endswith('.fa') or \
        file_lower.endswith('.fasta') or file_lower.endswith('.fna'):
            with open(self.fileFASTA, "r") as f:
                return self.ParseFASTA(f)
        else:
            print("Not in FASTA format.")

    def ParseFASTA(self, fileFASTA):
        '''Gets the sequence name and sequence from a FASTA formatted file'''
        fasta_list=[]
        for line in fileFASTA:
            if line[0] == '>':
                try:
                    fasta_list.append(current_dna)
            	#pass if an error comes up
                except UnboundLocalError:
                    #print "Inproper file format."
                    pass
                current_dna = [line.lstrip('>').rstrip('\n'),'']
            else:
                current_dna[1] += "".join(line.split())
        fasta_list.append(current_dna)
        '''Returns fasa as nested list, containing line identifier \
            and sequence'''
        return fasta_list



def clean_msa_names():
    tree = Phylo.read(mydir + 'data/tree/nmicrobiol201648-s8.txt', 'newick')
    read_fasta = classFASTA(mydir + 'data/align/nmicrobiol201648-s7.txt')
    OUT = open(mydir + 'data/align/nmicrobiol201648-s7_clean.txt','w+')
    tree_names = []

    for node in tree.get_terminals():
        tree_names.append(node.name)
    msa_names = []
    for x in read_fasta.readFASTA():
        header = x[0]
        for tup in to_replace:
            if tup[0] in header:
                header = header.replace(tup[0], tup[1])
        seq = x[1].replace('U', 'T')
        print>> OUT, '>' + header, '\n', seq, '\n'
        #msa_names.append(header)
    OUT.close()
    #tree_names_set = set(tree_names)
    #msa_names_set = set(msa_names)
    #intersect = tree_names_set & msa_names_set
    #tree_names_set_diff = tree_names_set - intersect
    #msa_names_set_diff = msa_names_set - intersect


def remove_bootstrap_values():

    #tree = Phylo.read(mydir + 'data/tree/nmicrobiol201648-s8.txt', 'newick')
    tree = open(mydir + 'data/tree/nmicrobiol201648-s8.txt','r')
    tree_out = open(mydir + 'data/tree/nmicrobiol201648-s8_clean.txt', 'w+')
    for line in tree:

        line_clean = re.sub(r'\[.*?\]', '', line)
        print>> tree_out, line_clean
    tree_out.close()
    #Phylo.write(tree, mydir + 'data/tree/nmicrobiol201648-s8_clean.txt', 'newick')
    #print tree

remove_bootstrap_values()
