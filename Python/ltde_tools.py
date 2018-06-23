import os, math
import pandas as pd

def get_path():
    return os.path.expanduser("~/GitHub/LTDE")

def clean_demography_df(df):
    df.loc[((df["strain"] == 'KBS0711W') & (df["rep"] == 1)), 'rep'] = 5
    df.loc[((df["strain"] == 'KBS0711W') & (df["rep"] == 2)), 'rep'] = 6
    df.loc[((df["strain"] == 'KBS0711W') & (df["rep"] == 3)), 'rep'] = 7
    df.loc[((df["strain"] == 'KBS0711W') & (df["rep"] == 4)), 'rep'] = 8
    df.loc[(df["strain"] == 'KBS0711W'), 'strain'] = 'KBS0711'
    return df


def get_strain_genus_dict():
    df_path = get_path() + '/data/traits/traits.txt'
    df = pd.read_csv(df_path, sep = '\t')
    genus_dict = pd.Series(df.Genus.values,index=df.Code).to_dict()
    genus_dict['KBS0727'] = 'Bradyrhizobium'
    return genus_dict

def weibull_mean(alpha, beta):
    return beta * math.gamma(1 + (1/alpha))

def weibull_half_life(alpha, beta):
    #return beta * math.gamma(1 + (1/alpha))
    return beta *  ((-math.log(0.5)) ** (1/alpha))

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


def get_mean_time_death():
    df_path = get_path() + '/data/demography/weibull_results.csv'
    df = pd.read_csv(df_path, sep = ',', index_col = 0)
    df['mean_days_death'] = df.apply(lambda row: weibull_mean(alpha = row['alpha'], beta = row['beta']), axis=1)
    df['sd_days_death'] = df.apply(lambda row: math.sqrt(weibull_variance(alpha = row['alpha'], beta = row['beta'])), axis=1)
    df['CI025_mean_days_death'] = df.apply(lambda row: weibull_CIs(mean = row['mean_days_death'], sd = row['sd_days_death'], n= row['N.obs'], lower = True ) , axis=1)
    df['CI975_mean_days_death'] = df.apply(lambda row: weibull_CIs(mean = row['mean_days_death'], sd = row['sd_days_death'], n= row['N.obs'], lower = False ) , axis=1)
    df['half_life'] = df.apply(lambda row: weibull_half_life(alpha = row['alpha'], beta = row['beta']), axis=1)
    return df



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
