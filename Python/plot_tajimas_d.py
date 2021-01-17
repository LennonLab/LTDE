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





df_taxa = pd.read_csv(lt.get_path() + '/data/breseq/tajimas_d_taxa.txt', sep = '\t')
df_taxa = df_taxa.sort_values(by=['tajimas_d'])
taxa_to_keep = df_taxa.Species.to_list()
df_taxa_samples = pd.read_csv(lt.get_path() + '/data/breseq/genetic_diversity.txt', sep = '\t', index_col=None)
fig = plt.figure()
for i, taxon in enumerate(taxa_to_keep):
    #print(taxon)
    #print(df_taxa_samples.loc[df_taxa_samples['Species'] == taxon])
    x_i = df_taxa_samples.loc[df_taxa_samples['Species'] == taxon].tajimas_d.values
    if len(x_i) < 3:
        taxa_to_keep.remove(taxon)
        continue
    x_i_mean = np.mean(x_i)
    x_i_sem = np.std(x_i)/ np.sqrt(len(x_i))
    taxon_color = lt.df_colors.loc[lt.df_colors['strain'] == taxon].Color.to_list()[0]

    jitter = np.random.normal(i, 0.04, size=len(x_i))
    plt.axvline(0, color = 'dimgrey', lw = 2, ls = '--', zorder=1)
    plt.errorbar(x_i_mean, i, xerr = x_i_sem*2, \
        fmt = 'o', alpha = 0.9, barsabove = True, marker = '.', \
        mfc = 'k', mec = 'k', c = 'k', zorder=3, ms=17)
    plt.scatter(x_i, jitter, c=taxon_color, marker = 'o', s = 80, \
        edgecolors='#244162', linewidth = 0, alpha = 1, zorder=2)

    p_value = df_taxa.loc[df_taxa['Species'] == taxon].p_BH.to_list()[0]
    if p_value < 0.05:
        plt.text(x_i_mean-0.1, i+0.15, r'$\ast$')


y1 = list(range(len(taxa_to_keep)))
latex_labels = [lt.latex_dict[x] for x in taxa_to_keep]

plt.xlabel("Tajima's D, " + r'$D_{T}$', fontsize = 12)

plt.yticks(y1, latex_labels, rotation=0)
fig.savefig(lt.get_path() + '/figs/tajimas_d.pdf', format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
