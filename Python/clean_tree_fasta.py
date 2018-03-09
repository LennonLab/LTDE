from __future__ import division
import os

mydir = os.path.expanduser("~/GitHub/LTDE/")


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



def rename_fasta():
    read_fasta = classFASTA(mydir + 'data/align/nmicrobiol201648-s7.txt')
    read_ltde_fasta = classFASTA(mydir + 'data/align/persistence_707_712_721.fasta')
    OUT = open(mydir + 'data/align/merged_align.txt','w+')
    for x in read_fasta.readFASTA():
        header = x[0].translate(None, '-.')
        header = header.replace(',', '_')
        header = header.replace(':', '_')
        header = header.replace('/', '_')
        header = header.replace("'", '_')
        seq = x[1].replace('U', 'T')
        seq = x[1].replace('-', '')
        #if header == 'Archaea_Diapherotrites_candidate_division_DUSEL3_archaeon_SCGC_AAA011E11_DUSEL_001_206':
        #    header = header[:8] + 'DPANN_' + header[8:]
        #if header == 'Bacteria_bacterium_JKG1':
        #    header = 'Bacteria_Chloroflexi_bacterium_JKG1'
        #if header == 'Bacteria_unclassified_Bacteria_Thermobaculum_Thermobaculum_terrenum_ATCC_BAA798':
        #    header = 'Bacteria_Chloroflexi_Thermobaculum_Thermobaculum_terrenum_ATCC_BAA798'
        #if 'Eukaryota' in header:
        #    continue
        print>> OUT, '>' + header, '\n', seq, '\n'

    for x in read_ltde_fasta.readFASTA():
        if x[0] == 'Methanosarcina':
            continue
        print>> OUT, '>' + x[0], '\n', x[1], '\n'

    OUT.close()

rename_fasta()
