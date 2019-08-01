from __future__ import division
import os, math, decimal
import numpy as np
import pandas as pd
from scipy.stats import poisson
from scipy.special import gammaln
from Bio import SeqIO

np.random.seed(123456789)

def get_path():
    return os.path.expanduser("~/GitHub/LTDE")


def get_strain_genus_dict():
    df_path = get_path() + '/data/traits/traits.txt'
    df = pd.read_csv(df_path, sep = '\t')
    genus_dict = pd.Series(df.Genus.values,index=df.Code).to_dict()
    genus_dict['KBS0727'] = 'Bradyrhizobium'
    return genus_dict

def get_genome_size_dict():
    genome_size_dict = {'KBS0703':4860540, 'ATCC13985':7045177, 'ATCC43928':6466929,
                        'KBS0701':6334343, 'KBS0702':3649547, 'KBS0705':4377150,
                        'KBS0706':8974938, 'KBS0707':6109947, 'KBS0710':6632121,
                        'KBS0711':6126717, 'KBS0712':7209490, 'KBS0713':4539645,
                        'KBS0714':2522703, 'KBS0715':3575377, 'KBS0721':5867608,
                        'KBS0722':4357994, 'KBS0724':7254422, 'KBS0725':7308410,
                        'KBS0727':7308448, 'KBS0801':7381868, 'KBS0802':6316887,
                        'KBS0812':4299846}

    return genome_size_dict

# NullMultiplicitySurvivalFunction class is modified from GitHub repo
# benjaminhgood/LTEE-metagenomic under GPL v2
class NullGeneMultiplicitySurvivalFunction(object):
    # Null multiplicity distribution for genes

    def __init__(self, Ls, ntot):
        self.ntot = ntot
        self.Ls = np.array(Ls)
        self.Lavg = self.Ls.mean()
        self.ps = self.Ls*1.0/self.Ls.sum()
        self.expected_ns = self.ntot*self.ps

    @classmethod
    def from_parallelism_statistics(cls, gene_parallelism_statistics):

        # calculate Ls
        Ls = []
        ntot = 0
        for gene_name in gene_parallelism_statistics.keys():
            Ls.append(gene_parallelism_statistics[gene_name]['length'])
            ntot += gene_parallelism_statistics[gene_name]['observed']

        return cls(Ls, ntot)

    def __call__(self, m):
        #lower_limits = np.ceil(m[:,None]*self.Ls[None,:]/self.Lavg)-1+0.1
        #return (poisson.sf(lower_limits, self.expected_ns[None,:])).sum(axis=1)
        lower_limits = np.ceil(m[:,None]*self.Ls[None,:]/self.Lavg)-2+0.1
        return (poisson.sf(lower_limits, self.expected_ns[None,:])*self.ps[None,:]).sum(axis=1)

# calculate_unnormalized_survival_from_vector function is modified from GitHub repo
# benjaminhgood/LTEE-metagenomic under GPL v2
def calculate_unnormalized_survival_from_vector(xs, min_x=None, max_x=None, min_p=1e-10):

    if min_x==None:
        min_x = xs.min()-1

    if max_x==None:
        max_x = xs.max()+1

    unique_xs = set(xs)
    unique_xs.add(min_x)
    unique_xs.add(max_x)

    xvalues = []
    num_observations = []

    for x in sorted(unique_xs):
        xvalues.append(x)
        num_observations.append( (xs>=x).sum() )

    # So that we can plot CDF, SF on log scale
    num_observations[0] -= min_p
    num_observations[1] -= min_p
    num_observations[-1] += min_p

    return np.array(xvalues), np.array(num_observations)

# calculate_G_scores function is modified from GitHub repo
# benjaminhgood/LTEE-metagenomic under GPL v2
def calculate_G_scores(gene_statistics, allowed_genes=None):
    # Calculates the G score for the whole gene, i.e.
    # n*g

    gene_g_scores = calculate_g_scores(gene_statistics,allowed_genes)

    gene_G_scores = {gene_name: gene_statistics[gene_name]['observed']*gene_g_scores[gene_name] for gene_name in gene_g_scores.keys()}

    return gene_G_scores

# calculate_total_parallelism function is modified from GitHub repo
# benjaminhgood/LTEE-metagenomic under GPL v2
def calculate_total_parallelism(gene_statistics, allowed_genes=None, num_bootstraps=10000):

    if allowed_genes==None:
        allowed_genes = gene_statistics.keys()

    Ls = []
    ns = []

    for gene_name in allowed_genes:

        Ls.append( gene_statistics[gene_name]['length'] )
        ns.append( gene_statistics[gene_name]['observed'] )


    Ls = np.array(Ls)
    ns = np.array(ns)

    Ltot = Ls.sum()
    ntot = ns.sum()
    ps = Ls*1.0/Ltot

    gs = ns*np.log(ns/(ntot*ps)+(ns==0))

    observed_G = gs.sum()/ns.sum()
    bootstrapped_Gs = []
    for bootstrap_idx in range(0,num_bootstraps):
        bootstrapped_ns = np.random.multinomial(ntot,ps)
        bootstrapped_gs = bootstrapped_ns*np.log(bootstrapped_ns/(ntot*ps)+(bootstrapped_ns==0))
        bootstrapped_G = bootstrapped_gs.sum()/bootstrapped_ns.sum()

        bootstrapped_Gs.append(bootstrapped_G)

    bootstrapped_Gs = np.array(bootstrapped_Gs)

    pvalue = ((bootstrapped_Gs>=observed_G).sum()+1.0)/(len(bootstrapped_Gs)+1.0)
    return observed_G, pvalue

# calculate_poisson_log_survival function is modified from GitHub repo
# benjaminhgood/LTEE-metagenomic under GPL v2
def calculate_poisson_log_survival(ns, expected_ns):

    survivals = poisson.sf(ns-0.1, expected_ns)

    logsurvivals = np.zeros_like(survivals)
    logsurvivals[survivals>1e-20] = -np.log(survivals[survivals>1e-20])
    logsurvivals[survivals<=1e-20] = (-ns*np.log(ns/expected_ns+(ns==0))+ns-expected_ns)[survivals<=1e-20]

    return logsurvivals

# calculate_parallelism_logpvalues function is modified from GitHub repo
# benjaminhgood/LTEE-metagenomic under GPL v2
def calculate_parallelism_logpvalues(gene_statistics):

    gene_names = []
    Ls = []
    ns = []
    expected_ns = []

    for gene_name in gene_statistics.keys():
        gene_names.append(gene_name)
        ns.append(gene_statistics[gene_name]['observed'])
        expected_ns.append(gene_statistics[gene_name]['expected'])

    ns = np.array(ns)
    expected_ns = np.array(expected_ns)

    #print(list(zip(ns, expected_ns)))

    logpvalues = calculate_poisson_log_survival(ns, expected_ns)

    return {gene_name: logp for gene_name, logp in zip(gene_names, logpvalues)}


# NullGeneLogpSurvivalFunction class is modified from GitHub repo
# benjaminhgood/LTEE-metagenomic under GPL v2
class NullGeneLogpSurvivalFunction(object):
    # Null distribution of -log p for each gene

    def __init__(self, Ls, ntot,nmin=0):
        self.ntot = ntot
        self.Ls = np.array(Ls)*1.0
        self.Lavg = self.Ls.mean()
        self.ps = self.Ls/self.Ls.sum()
        self.expected_ns = self.ntot*self.ps
        self.nmin = nmin

    @classmethod
    def from_parallelism_statistics(cls, gene_parallelism_statistics,nmin=0):

        # calculate Ls
        Ls = []
        ntot = 0
        for gene_name in gene_parallelism_statistics.keys():
            Ls.append(gene_parallelism_statistics[gene_name]['length'])
            ntot += gene_parallelism_statistics[gene_name]['observed']

        return cls(Ls, ntot, nmin)

    def __call__(self, mlogps):

        # Do sum by hand
        ns = np.arange(0,400)*1.0

        logpvalues = calculate_poisson_log_survival(ns[None,:], self.expected_ns[:,None])

        logprobabilities = ns[None,:]*np.log(self.expected_ns)[:,None]-gammaln(ns+1)[None,:]-self.expected_ns[:,None]
        probabilities = np.exp(logprobabilities)
        survivals = np.array([ ((logpvalues>=mlogp)*(ns[None,:]>=self.nmin)*probabilities).sum() for mlogp in mlogps])
        return survivals



# calculate_synonymous_nonsynonymous_target_sizes function is modified from GitHub repo
# benjaminhgood/LTEE-metagenomic under GPL v2
def calculate_synonymous_nonsynonymous_target_sizes(taxon):
    codon_table = { 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'CGT': 'R',
                    'CGC': 'R', 'CGA':'R', 'CGG':'R', 'AGA':'R', 'AGG':'R',
                    'AAT':'N', 'AAC':'N', 'GAT':'D', 'GAC':'D', 'TGT':'C',
                    'TGC':'D', 'CAA':'Q', 'CAG':'Q', 'GAA':'E', 'GAG':'E',
                    'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G', 'CAT':'H',
                    'CAC':'H', 'ATT':'I', 'ATC':'I', 'ATA':'I', 'TTA':'L',
                    'TTG':'L', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
                    'AAA':'K', 'AAG':'K', 'ATG':'M', 'TTT':'F', 'TTC':'F',
                    'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'TCT':'S',
                    'TCC':'S', 'TCA':'S', 'TCG':'S', 'AGT':'S', 'AGC':'S',
                    'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T', 'TGG':'W',
                    'TAT':'Y', 'TAC':'Y', 'GTT':'V', 'GTC':'V', 'GTA':'V',
                    'GTG':'V', 'TAA':'!', 'TGA':'!', 'TAG':'!' }

    # calculate number of synonymous opportunities for each codon
    codon_synonymous_opportunity_table = {}
    for codon in codon_table.keys():
        codon_synonymous_opportunity_table[codon] = {}
        for i in range(0,3):
            codon_synonymous_opportunity_table[codon][i] = -1 # since G->G is by definition synonymous, but we don't want to count it
            codon_list = list(codon)
            for base in ['A','C','T','G']:
                codon_list[i]=base
                new_codon = "".join(codon_list)
                if codon_table[codon]==codon_table[new_codon]:
                    # synonymous!
                    codon_synonymous_opportunity_table[codon][i]+=1

    codon_synonymous_substitution_table = {}
    codon_nonsynonymous_substitution_table = {}
    for codon in codon_table.keys():
        codon_synonymous_substitution_table[codon] = [[],[],[]]
        codon_nonsynonymous_substitution_table[codon] = [[],[],[]]

        for i in range(0,3):
            reference_base = codon[i]

            codon_list = list(codon)
            for derived_base in ['A','C','T','G']:
                if derived_base==reference_base:
                    continue
                substitution = '%s->%s' % (reference_base, derived_base)
                codon_list[i]=derived_base
                new_codon = "".join(codon_list)
                if codon_table[codon]==codon_table[new_codon]:
                    # synonymous!
                    codon_synonymous_substitution_table[codon][i].append(substitution)
                else:
                    codon_nonsynonymous_substitution_table[codon][i].append(substitution)
    bases = set(['A','C','T','G'])
    substitutions = []
    for b1 in bases:
        for b2 in bases:
            if b2==b1:
                continue
            substitutions.append( '%s->%s' % (b1,b2) )
    substitution_specific_synonymous_sites = {substitution: 0 for substitution in substitutions}
    substitution_specific_nonsynonymous_sites = {substitution: 0 for substitution in substitutions}
    effective_gene_synonymous_sites = {}
    effective_gene_nonsynonymous_sites = {}
    gene_length_map = {}
    genome_path = get_path() + '/data/genomes/genomes_ncbi/' + taxon
    for subdir, dirs, files in os.walk(genome_path):
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

                            if gene_name not in effective_gene_synonymous_sites:
                                effective_gene_synonymous_sites[gene_name]=0
                                effective_gene_nonsynonymous_sites[gene_name]=0
                            aa_str = str(feature.qualifiers['translation'][0])
                            nuc_str = str(feature.location.extract(record).seq[:-3])
                            gene_length_map[gene_name] = len(nuc_str)
                            for position in range(len(nuc_str)):
                                codon_start = int(position/3)*3
                                codon = nuc_str[codon_start:codon_start+3]
                                if len(codon) <3:
                                    continue
                                position_in_codon = position%3

                                effective_gene_synonymous_sites[gene_name] += codon_synonymous_opportunity_table[codon][position_in_codon]/3.0
                                effective_gene_nonsynonymous_sites[gene_name] += 1-codon_synonymous_opportunity_table[codon][position_in_codon]/3.0

                                for substitution in codon_synonymous_substitution_table[codon][position_in_codon]:
                                    substitution_specific_synonymous_sites[substitution] += 1

                                for substitution in codon_nonsynonymous_substitution_table[codon][position_in_codon]:
                                    substitution_specific_nonsynonymous_sites[substitution] += 1

    substitution_specific_synonymous_fraction = {substitution: substitution_specific_synonymous_sites[substitution]*1.0/(substitution_specific_synonymous_sites[substitution]+substitution_specific_nonsynonymous_sites[substitution]) for substitution in substitution_specific_synonymous_sites.keys()}
    effective_gene_lengths = {gene_name: gene_length_map[gene_name]-effective_gene_synonymous_sites[gene_name] for gene_name in gene_length_map.keys()}
    effective_gene_lengths_synonymous = sum([effective_gene_synonymous_sites[gene_name] for gene_name in gene_length_map.keys()])
    effective_gene_lengths_nonsynonymous = sum([effective_gene_nonsynonymous_sites[gene_name] for gene_name in gene_length_map.keys()])
    effective_gene_lengths_noncoding = get_genome_size_dict()[taxon] - effective_gene_lengths_synonymous-effective_gene_lengths_nonsynonymous

    return effective_gene_lengths, effective_gene_lengths_synonymous, effective_gene_lengths_nonsynonymous, effective_gene_lengths_noncoding





def get_final_pop_size(population):
    taxon = population.split('-')[0]
    rep = rename_rep()[population.split('-')[1]]
    df = pd.read_csv(get_path() + '/data/demography/longtermdormancy_20190528_nocomments.csv', sep = ',')
    df_taxon_rep = df[(df['Strain']==taxon) & (df['Rep']==rep)]
    df_taxon_rep['Firstread_date'] = pd.to_datetime(df_taxon_rep['Firstread_date'])
    last_sample = df_taxon_rep.loc[df_taxon_rep['Firstread_date'] == df_taxon_rep['Firstread_date'].max()]
    return int((int(last_sample.Colonies)+1) * (1000/float(last_sample.Inoculum)) * (10**float(last_sample.Dilution)))



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



def get_pooled_pi(nested_list, size):
    pi_sum = 0
    for i in nested_list:
        pi_sum += (2 / (i[1] * (i[1]-1)))  * (i[3] * (i[1]-i[3]))
    return pi_sum / size


def get_pi(nested_list, n_c, size):
    pi_sum = 0
    for i in nested_list:
        pi_sum += 2*i[0]*(1-i[0])
    return (n_c/ (n_c-1)) * (pi_sum / size)

def get_a_1(n_c):
    return sum( [1/k for k in range(1, n_c)] )

def get_theta(nested_list, n_c, size):
    return len(nested_list[0]) / (get_a_1(n_c) * size)



def get_TD(freq_list, pi, theta, n_c):
    S = len(freq_list[0])
    a_1 = get_a_1(n_c)
    b_1 = (n_c+1) / (3*(n_c-1))
    c_1 = b_1 - (1/a_1)
    e_1 = c_1/a_1

    a_2 = sum( [1/(k**2) for k in range(1, n_c)] )
    b_2 = 2*((n_c**2)+n_c+3) / (9*n_c*(n_c-1))
    c_2 = b_2 - ((n_c+2)/(a_1*n_c)) + (a_2/(a_1**2))
    e_2 = c_2 /((a_1 **2) + a_2)

    return (pi - theta) / math.sqrt( (e_1*S)+(e_2*S*(S-1)) )
