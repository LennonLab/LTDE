import os

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

def split_by_n( seq, n ):
    """A generator to divide a sequence into chunks of n units."""
    while seq:
        yield seq[:n]
        seq = seq[n:]

def rename_fasta(IN, OUT):
    read_FASTA = classFASTA(IN)
    OUT_FASTA = open(OUT,'w+')
    for x in read_FASTA.readFASTA():
        #x[0].strip #### include this command
        print>> OUT_FASTA, '>' + x[0].strip().replace('|', '_')
        split_seq = split_by_n(x[1], 60)
        for seq in split_seq:
            print>> OUT_FASTA, seq
        # '\n', x[1], '\n'
    OUT_FASTA.close()

#def rename_fasta_all():

IN = '/N/dc2/projects/muri2/Task2/LTDE/data/reference_genomes/2016_KBSGenomes_Annotate/ATCC13985/G-Chr1.fna'
OUT = '/N/dc2/projects/muri2/Task2/LTDE/data/reference_genomes/genomes_rename_fna/ATCC13985.fna'

#IN = '/Users/WRShoemaker/GitHub/LTDE/data/genomes/2016_KBSGenomes_Annotate/ATCC13985/G-Chr1.fna'
#OUT = '/Users/WRShoemaker/GitHub/LTDE/data/genomes/genomes_rename_fna/ATCC13985.fna'


rename_fasta(IN, OUT)

#'gnl|I|G_contig000058' to 'gnl_I_G_contig000058'
