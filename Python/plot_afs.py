from __future__ import division

import glob, math, os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import ltde_tools as lt
from scipy import stats
from scipy.stats import t
from scipy.integrate import odeint
from decimal import Decimal
import _pickle as pickle
from matplotlib.ticker import FormatStrFormatter
from matplotlib.lines import Line2D

import matplotlib.lines as mlines

import matplotlib.ticker
import datetime as dt

#from sklearn.model_selection import GridSearchCV
#from sklearn.neighbors import KernelDensity
import statsmodels.stats.multitest as mt
import statsmodels.formula.api as smf

from Bio import SeqIO

from statsmodels.base.model import GenericLikelihoodModel

# only plot taxa w/ significant g scores and at least 100 mutations


n_bins=20

afs_taxa_reps = [x.split('.')[0] for x in os.listdir(lt.get_path() + '/data/breseq/allele_freq_spec/') if x.endswith(".txt")]
afs_taxa = list(set([x.split('.')[0].split('-')[0] for x in os.listdir(lt.get_path() + '/data/breseq/allele_freq_spec/') if x.endswith(".txt")]))
_nsre = re.compile('([0-9]+)')
def natural_sort_key(s):
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split(_nsre, s)]
afs_taxa.sort(key=natural_sort_key)
taxa_to_analyze = []
for taxon in afs_taxa:
    taxon_reps = [x for x in afs_taxa_reps if taxon in x]
    if len(taxon_reps) >= 3:
        taxa_to_analyze.append(taxon)
#print(str(len(taxa_to_analyze)) +" taxa to analyze")
fig = plt.figure()
#taxa_to_analyze = taxa_to_analyze[:2]
fig.subplots_adjust(hspace=0.35, wspace=0.35)

for i, taxon in enumerate(lt.taxa_to_plot):
    taxon_color = lt.df_colors.loc[lt.df_colors['strain'] == taxon].Color.to_list()[0]
    taxon_reps = [x.split('-')[1] for x in afs_taxa_reps if taxon in x]
    taxon_reps.sort()
    ax = fig.add_subplot(3, 3, i+1)
    taxon_means = []
    taxons_max = []
    freqs_all = []
    for taxon_rep in taxon_reps:
        df = pd.read_csv(lt.get_path() + '/data/breseq/allele_freq_spec/' + taxon + '-'+taxon_rep +'.txt', sep = '\t')
        freqs = df.frequency.tolist()
        freqs_all.extend(freqs)
        #grid_ = GridSearchCV(KernelDensity(),
        #                {'bandwidth': np.linspace(0.001, 1, 50)},
        #                cv=2) # 20-fold cross-validation
        #grid_.fit(freqs[:, None])
        #x_grid_ = np.linspace(-2, 0, 1000)
        #x_grid_ = np.linspace(0, 1, 1000)
        #print(grid_.best_params_)
        #kde_ = grid_.best_estimator_
        #pdf_ = np.exp(kde_.score_samples(x_grid_[:, None]))
        #pdf_ = kde_.score_samples(x_grid_[:, None])
        #pdf_ = [x / sum(pdf_) for x in pdf_]


        #ax.hist(freqs_all, bins= np.logspace(np.log10(0.1),np.log10(1.0), 50), alpha=0.8,  color = taxon_color, weights=np.zeros_like(freqs) + 1. / len(freqs))
        #ax.plot(x_grid_, pdf_, alpha=0.8, lw = 2, color = taxon_color) #, marker='o')
        #ax.plot(10**x_grid_, pdf_, alpha=0.8, lw = 2, color = taxon_color) #, marker='o')
        taxon_means.append(np.mean(freqs))
        taxons_max.append(max(freqs))
        #ax.set_xscale('log')

    freqs_all = np.asarray(freqs_all)
    #hist, bins, _ = plt.hist(freqs_all, bins=40)
    #hist, bins = np.histogram(freqs_all, bins=n_bins)
    #logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
    # histogram on log scale.
    # Use non-equal bin sizes, such that they look equal on log scale.
    #logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))

    ax.hist(freqs_all, bins=40, alpha=1,  color = taxon_color, weights=np.zeros_like(freqs_all) + 1. / len(freqs_all))

    ax.title.set_text(lt.latex_dict[taxon])
    ax.title.set_fontsize(7)

    #ax.axvline(x=np.mean(taxon_means), color='black', linestyle='--', lw = 1.5, alpha=0.8)
    #ax.axvline(x=np.mean(taxons_max), color='red', linestyle='--', lw = 1.5, alpha=0.8)

    ax.axvline(x=np.mean(taxon_means), color='black', linestyle='--', lw = 1.5, alpha=0.8)
    ax.axvline(x=np.mean(taxons_max), color='red', linestyle='--', lw = 1.5, alpha=0.8)

    ax.tick_params(axis='both', which='major', labelsize=7)
    ax.tick_params(axis='both', which='minor', labelsize=5)

    ax.set_xlim([min(freqs_all), 0.65])

    #ax.set_xscale('log', basex=10)

fig.subplots_adjust(wspace=0.3, hspace=0.4)
fig.text(0.5, 0.02, 'Mutation frequency', ha='center', fontsize=16)
fig.text(0.02, 0.5, 'Frequency', va='center', rotation='vertical', fontsize=16)

fig.savefig(lt.get_path() + '/figs/afs.pdf', format = 'pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
