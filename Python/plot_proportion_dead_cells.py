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




df = pd.read_csv(lt.get_path() + '/data/staining.all.new.txt', sep = '\t')
df = df[df.strain != "KBS0725"]
taxa = list(set(df.strain.to_list()))
to_remove = ['KBS0711W', 'KBS0727', 'KBS0714', 'KBS0701']
taxa = [ x for x in taxa if x not in to_remove]
df_anc = df.loc[df['hist'] == 'anc']
df_der = df.loc[df['hist'] == 'der']
fig = plt.figure()
plt.axvline(1, color = 'dimgrey', lw = 2, ls = '--', zorder=1)
# first sort by mean
mean_list = []
for taxon in taxa:
    df_der_dead = df_der.loc[df_der['strain'] == taxon].dead.values
    delta_dead = df_der_dead - df_anc.loc[df_anc['strain'] == taxon].dead.values[0]
    delta_dead_mean = np.mean(delta_dead)
    mean_list.append((taxon, delta_dead_mean))

mean_list_sorted = sorted(mean_list, key=lambda tup: tup[1])
taxa_sorted = [x[0] for x in mean_list_sorted][::-1]
p_values = []
for i, taxon in enumerate(taxa_sorted):
    taxon_color = lt.df_colors.loc[lt.df_colors['strain'] == taxon].Color.to_list()[0]
    df_der_dead = df_der.loc[df_der['strain'] == taxon].dead.values
    delta_dead = df_der_dead - df_anc.loc[df_anc['strain'] == taxon].dead.values[0]

    delta_dead_mean = np.mean(delta_dead)
    delta_dead_sem = np.std(delta_dead)/ np.sqrt(len(delta_dead))

    jitter = np.random.normal(i, 0.04, size=len(delta_dead))

    plt.errorbar(delta_dead_mean, i, xerr = delta_dead_sem*2, \
        fmt = 'o', alpha = 0.9, barsabove = True, marker = '.', \
        mfc = 'k', mec = 'k', c = 'k', zorder=3, ms=17)
    plt.scatter(delta_dead, jitter, c=taxon_color, marker = 'o', s = 80, \
        edgecolors='#244162', linewidth = 0, alpha = 1, zorder=2)


y1 = list(range(len(taxa_sorted)))
latex_labels = [lt.latex_dict[x] for x in taxa]


plt.xlabel('Change in proportion of dead cells', fontsize = 12)
plt.xlim(-1,1.2)
plt.yticks(y1, latex_labels, rotation=0)
fig.savefig(lt.get_path() + '/figs/prop_dead.pdf', format= 'pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
