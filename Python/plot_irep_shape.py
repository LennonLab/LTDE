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



df_irep = pd.read_csv(lt.get_path() + '/data/iRep_clean.txt', sep = '\t')
df_irep = df_irep.rename(columns={'Species': 'strain'})
df_weib = pd.read_csv(lt.get_path() + '/data/demography/weibull_results_clean.csv', sep = ',')
df_merged = df_weib.merge(df_irep, on=['strain','rep'])
taxa = list(set(df_merged.strain.to_list()))


df_merged['alpha_log10'] = np.log10(df_merged.alpha)

mf = smf.mixedlm("alpha_log10 ~ iRep", df_merged, groups=df_merged["strain"])
mf_fit = mf.fit()
print(mf_fit.summary())

irep_mean_list = []
shape_mean_list = []

fig = plt.figure()
plt.axhline(1, color = 'dimgrey', lw = 2, ls = '--', zorder=1)
for i, taxon in enumerate(taxa):

    taxon_color = lt.df_colors.loc[lt.df_colors['strain'] == taxon].Color.to_list()[0]
    irep = df_merged.loc[df_merged['strain'] == taxon].iRep.values
    shape = df_merged.loc[df_merged['strain'] == taxon].alpha.values

    irep_mean = np.mean(irep)
    irep_sem = np.std(irep)/ np.sqrt(len(irep))
    irep_mean_list.append(irep_mean)

    shape_mean = 10**np.mean(np.log10(shape))
    shape_sem = (np.std(np.log10(shape))/ np.sqrt(len(irep)))
    shape_mean_list.append(np.mean(np.log10(shape)))

    plt.errorbar(irep_mean, shape_mean, xerr = irep_sem*2, yerr = shape_sem*2, \
        fmt = 'o', alpha = 0.9, barsabove = True, marker = '.', \
        mfc = 'k', mec = 'k', c = 'k', zorder=2, ms=17)
    plt.scatter(irep_mean, shape_mean, c=taxon_color, marker = 'o', s = 93, \
        edgecolors='none', linewidth = 0, alpha = 1, zorder=3)


    #plt.errorbar(irep_mean, shape_mean, xerr= alpha_sem, yerr = dn_ds_sem*2, \
    #    fmt = 'o', alpha = 0.9, barsabove = True, marker = '.', \
    #    mfc = 'none', mec = 'none', ecolor = 'k', zorder=1, ms=15)
    #plt.scatter(alpha_mean, dn_ds_mean, c=taxon_color, s=17, zorder=2)


plt.ylim(0.03, 1.1)
plt.yscale('log',basey=10)
plt.xlabel('Index of replication (iRep)', fontsize = 12)
plt.ylabel('Shape paramter, ' r'$k$', fontsize = 12)
fig.savefig(lt.get_path() + '/figs/irep_shape.pdf', format = 'pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()

slope, intercept, r_value, p_value, std_err = stats.linregress(irep_mean_list, shape_mean_list)


print("Mean iRep across taxa = " + str(round(np.mean(irep_mean_list), 4)))

print("D.F. = " + str(len(irep_mean_list)-2))
print("t = " + str(round((slope-0)/std_err, 4)))
print("r^2 = " + str(round(r_value**2, 4)))
print("p-value = " + str(round(p_value, 4)))
