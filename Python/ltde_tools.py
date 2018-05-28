import os, math

def get_path():
    return os.path.expanduser("~/GitHub/LTDE")

def weibull_mean(alpha, beta):
    return beta * math.gamma(1 + (1/alpha))

def weibull_variance(alpha, beta):
    return (beta ** 2) * (math.gamma(1 + (2/alpha)) -  (math.gamma(1 + (1/alpha)) ** 2))

def weibull_CIs(mean, sd, n, lower = True, pooled = False):
    # alpha of 0.05, so z_(alpha/2) = 1.96
    if pooled == True:
        # pass list of sample sizes
        n_pooled = math.sqrt( sum(1 / n) )
        se_g1 = (sd * n_pooled) / mean
    else:
        se_g1 = (sd / math.sqrt(n)) / mean
    if lower == True:
        CI = math.exp( math.log(mean) - (1.96 * se_g1))
    else:
        CI = math.exp( math.log(mean) + (1.96 * se_g1))
    return CI


class classFASTA:

    def __init__(self, fileFASTA):
        self.fileFASTA = fileFASTA

    def readFASTA(self):
        '''Checks for fasta by file extension'''
        file_lower = self.fileFASTA.lower()
        '''Check for three most common fasta file extensions'''
        if file_lower.endswith('.txt') or file_lower.endswith('.fa') or \
        file_lower.endswith('.fasta') or file_lower.endswith('.fna') or \
        file_lower.endswith('.faa') or file_lower.endswith('.ffn'):
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


def strain_list():
    strains = ['ATCC13985', 'ATCC43928', 'KBS0701', 'KBS0702', 'KBS0703', 'KBS0705', \
                'KBS0706', 'KBS0707', 'KBS0710', 'KBS0711', 'KBS0712', 'KBS0713', 'KBS0714', \
                'KBS0715', 'KBS0721', 'KBS0722', 'KBS0724', 'KBS0725', 'KBS0727', 'KBS0801', \
                'KBS0802', 'KBS0812']
    return strains

def rename_rep():
    return {'A':1, 'B':2, 'C':3, 'C1':3, 'D':4, 'E':5, 'F':6, 'K':10, 'L':11}

def split_by_n( seq, n ):
    """A generator to divide a sequence into chunks of n units."""
    while seq:
        yield seq[:n]
        seq = seq[n:]
