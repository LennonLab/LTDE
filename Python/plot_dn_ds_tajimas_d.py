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





df_taxa = pd.read_csv(lt.get_path() + '/data/breseq/dN_dS_taxa.txt', sep = '\t')
df_taxa = df_taxa.sort_values(by=['dN_dS_total'])
taxa_to_keep = df_taxa.Species.to_list()
df_diversity = pd.read_csv(lt.get_path() + '/data/breseq/genetic_diversity.txt', sep = '\t', index_col=None)
df_diversity['mean_birth_per_death_log10'] = np.log10(df_diversity['mean_birth_per_death'].values)

df_demographic = pd.read_csv(lt.get_path() + '/data/demography/weibull_results_clean.csv', sep = ',', index_col=None)
df_demographic.rename(columns={'strain':'Species'}, inplace=True)

df_merged = pd.merge(df_diversity, df_demographic,  how='left', on=['Species','rep'])#, right_on = ['Species','c2'])
df_merged['N_0_beta_log10'] = np.log10(df_merged['N_0_beta'].values)
df_merged['alpha_log10'] = np.log10(df_merged['alpha'].values)
df_merged['N_delta_log10'] = np.log10(df_merged['N_0'].values - df_merged['N_final'].values)
df_merged['mean_binary_divisions_log10'] = np.log10(df_merged['mean_binary_divisions'].values)


print("N_delta log10 vs dNdS ")
N_delta_dnds = smf.mixedlm("dn_ds_total ~ N_delta_log10", df_merged, groups=df_merged["Species"])
N_delta_dnds_fit = N_delta_dnds.fit()
print(N_delta_dnds_fit.summary())

print()

print("N_0_beta_log10 vs dNdS ")
N_0_beta_dnds = smf.mixedlm("dn_ds_total ~ N_0_beta_log10", df_merged, groups=df_merged["Species"])
N_0_beta_dnds_fit = N_0_beta_dnds.fit()
print(N_0_beta_dnds_fit.summary())

print()

print("alpha vs dNdS ")
alpha_dnds = smf.mixedlm("dn_ds_total ~ alpha", df_merged, groups=df_merged["Species"])
alpha_dnds_fit = alpha_dnds.fit()
print(alpha_dnds_fit.summary())


print()

print("mean_binary_divisions_log10 vs dNdS ")
b_dnds = smf.mixedlm("dn_ds_total ~ mean_binary_divisions_log10", df_merged, groups=df_merged["Species"])
b_dnds_fit = b_dnds.fit()
print(b_dnds_fit.summary())


print()


print("N_delta log10 vs Tajimas D  ")
N_delta_TD = smf.mixedlm("tajimas_d ~ N_delta_log10", df_merged, groups=df_merged["Species"])
N_delta_TD_fit = N_delta_TD.fit()
print(N_delta_TD_fit.summary())

print()

print("N_0_beta_log10 vs Tajimas D  ")
N_0_beta_TD = smf.mixedlm("tajimas_d ~ N_0_beta_log10", df_merged, groups=df_merged["Species"])
N_0_beta_TD_fit = N_0_beta_TD.fit()
print(N_0_beta_TD_fit.summary())

print()

print("alpha vs Tajimas D  ")
alpha_TD = smf.mixedlm("tajimas_d ~ alpha", df_merged, groups=df_merged["Species"])
alpha_TD_fit = alpha_TD.fit()
print(alpha_TD_fit.summary())


print()

print("mean_binary_divisions_log10 vs Tajimas D ")
b_TD = smf.mixedlm("tajimas_d ~ mean_binary_divisions_log10", df_merged, groups=df_merged["Species"])
b_TD_fit = b_TD.fit()
print(b_TD_fit.summary())



#print(t_N_delta_fit.params)

#covb = t_N_delta_fit.cov_params()
#prediction_var = t_N_delta_fit.mse_resid + (df_merged.mean_binary_divisions_log10.values * np.dot(covb,df_merged.mean_binary_divisions_log10.values.T).T).sum(1)

#print(t_N_delta.predict(t_N_delta_fit.params))

#print(t_N_delta.get_prediction(df_merged.mean_binary_divisions_log10.values))
#fixed effects variance
#print(t_N_delta_fit.params)
#f_var = np.var(t_N_delta.predict(t_N_delta_fit.params))


fig = plt.figure(figsize = (6, 3))
fig.tight_layout(pad = 2.8)

ax1 = plt.subplot2grid((2, 4), (0, 0), colspan=1)
ax2 = plt.subplot2grid((2, 4), (0, 1), colspan=1)
ax3 = plt.subplot2grid((2, 4), (0, 2), colspan=1)
ax4 = plt.subplot2grid((2, 4), (0, 3), colspan=1)

ax5 = plt.subplot2grid((2, 4), (1, 0), colspan=1)
ax6 = plt.subplot2grid((2, 4), (1, 1), colspan=1)
ax7 = plt.subplot2grid((2, 4), (1, 2), colspan=1)
ax8 = plt.subplot2grid((2, 4), (1, 3), colspan=1)

for i, taxon in enumerate(taxa_to_keep):
    df_taxon_samples = df_merged.loc[df_merged['Species'] == taxon]

    dn_ds = df_merged.loc[df_merged['Species'] == taxon].dn_ds_total.values
    dn_ds_mean = np.mean(dn_ds)
    dn_ds_sem = np.std(dn_ds)/ np.sqrt(len(dn_ds))

    t_d = df_merged.loc[df_merged['Species'] == taxon].tajimas_d.values
    t_d_mean = np.mean(t_d)
    t_d_sem = np.std(t_d)/ np.sqrt(len(t_d))

    N_delta_log10 = df_merged.loc[df_merged['Species'] == taxon].N_delta_log10.values
    N_delta_mean = 10 ** np.mean(N_delta_log10)
    N_delta_sem = 10 ** (2*np.std(N_delta_log10)/ np.sqrt(len(N_delta_log10)))

    alpha = df_merged.loc[df_merged['Species'] == taxon].alpha.values
    alpha_mean = np.mean(alpha)
    alpha_sem = (2*np.std(alpha)/ np.sqrt(len(alpha)))

    N_0_beta_log10 = df_merged.loc[df_merged['Species'] == taxon].N_0_beta_log10.values
    N_0_beta_mean = 10 ** np.mean(N_0_beta_log10)
    N_0_beta_sem = 10 ** (2*np.std(N_0_beta_log10)/ np.sqrt(len(alpha)))

    N_births = df_merged.loc[df_merged['Species'] == taxon].mean_binary_divisions.values
    N_births_mean = 10 ** np.mean(np.log10(N_births))
    N_births_sem = 10 ** (2*np.std(np.log10(N_births))/ np.sqrt(len(alpha)))

    taxon_color = lt.df_colors.loc[lt.df_colors['strain'] == taxon].Color.to_list()[0]

    ax1.errorbar(N_delta_mean, dn_ds_mean, xerr= N_delta_sem, yerr = dn_ds_sem*2, \
        fmt = 'o', alpha = 0.9, barsabove = True, marker = '.', \
        mfc = 'none', mec = 'none', ecolor = 'k', zorder=1, ms=15)
    ax1.scatter(N_delta_mean, dn_ds_mean, c=taxon_color, s=12, zorder=2)


    ax2.errorbar(N_0_beta_mean, dn_ds_mean, xerr= N_0_beta_sem, yerr = dn_ds_sem*2, \
        fmt = 'o', alpha = 0.9, barsabove = True, marker = '.', \
        mfc = 'none', mec = 'none', ecolor = 'k', zorder=1, ms=15)
    ax2.scatter(N_0_beta_mean, dn_ds_mean, c=taxon_color, s=12, zorder=2)


    ax3.errorbar(alpha_mean, dn_ds_mean, xerr= alpha_sem, yerr = dn_ds_sem*2, \
        fmt = 'o', alpha = 0.9, barsabove = True, marker = '.', \
        mfc = 'none', mec = 'none', ecolor = 'k', zorder=1, ms=15)
    ax3.scatter(alpha_mean, dn_ds_mean, c=taxon_color, s=12, zorder=2)


    ax4.errorbar(N_births_mean, dn_ds_mean, xerr= N_births_sem, yerr = dn_ds_sem*2, \
        fmt = 'o', alpha = 0.9, barsabove = True, marker = '.', \
        mfc = 'none', mec = 'none', ecolor = 'k', zorder=1, ms=15)
    ax4.scatter(N_births_mean, dn_ds_mean, c=taxon_color, s=12, zorder=2)



    ax5.errorbar(N_delta_mean, t_d_mean, xerr= N_delta_sem, yerr = t_d_sem*2, \
        fmt = 'o', alpha = 0.9, barsabove = True, marker = '.', \
        mfc = 'none', mec = 'none', ecolor = 'k', zorder=1, ms=15)
    ax5.scatter(N_delta_mean, t_d_mean, c=taxon_color, s=12, zorder=2)


    ax6.errorbar(N_0_beta_mean, t_d_mean, xerr= N_0_beta_sem, yerr = t_d_sem*2, \
        fmt = 'o', alpha = 0.9, barsabove = True, marker = '.', \
        mfc = 'none', mec = 'none', ecolor = 'k', zorder=1, ms=15)
    ax6.scatter(N_0_beta_mean, t_d_mean, c=taxon_color, s=12, zorder=2)


    ax7.errorbar(alpha_mean, t_d_mean, xerr= alpha_sem, yerr = t_d_sem*2, \
        fmt = 'o', alpha = 0.9, barsabove = True, marker = '.', \
        mfc = 'none', mec = 'none', ecolor = 'k', zorder=1, ms=15)
    ax7.scatter(alpha_mean, t_d_mean, c=taxon_color, s=12, zorder=2)


    ax8.errorbar(N_births_mean, t_d_mean, xerr= N_births_sem, yerr = t_d_sem*2, \
        fmt = 'o', alpha = 0.9, barsabove = True, marker = '.', \
        mfc = 'none', mec = 'none', ecolor = 'k', zorder=1, ms=15)
    ax8.scatter(N_births_mean, t_d_mean, c=taxon_color, s=12, zorder=2)



#p<0.0001
x_ax8= np.linspace(4, 5.6, num=500)
y_pred_ax8 = b_TD_fit.params[0] + x_ax8 * b_TD_fit.params[1]
ax8.plot(10**x_ax8, y_pred_ax8, c='k')

ax8.text(0.07, 0.9, r'$\beta_{1}=$' + str(round(b_TD_fit.params[1], 2)), fontsize=5, transform=ax8.transAxes)
ax8.text(0.07, 0.8, r'$p\ll0.001$', fontsize=5, transform=ax8.transAxes)



ax7.text(0.07, 0.9, r'$\beta_{1}=$' + str(round(alpha_TD_fit.params[1], 2)), fontsize=5, transform=ax7.transAxes)
ax7.text(0.07, 0.8, r'$p=$' + str(round(alpha_TD_fit.pvalues[1], 3)), fontsize=5, transform=ax7.transAxes)

x_ax7 = np.linspace(0.2, 0.65, num=500)
y_pred_ax7 = alpha_TD_fit.params[0] + x_ax7 * alpha_TD_fit.params[1]
ax7.plot(x_ax7, y_pred_ax7, c='k')



ax6.text(0.6, 0.9, r'$\beta_{1}=$' + str(round(N_0_beta_TD_fit.params[1], 2)), fontsize=5, transform=ax6.transAxes)
ax6.text(0.6, 0.8, r'$p=$' + str(round(N_0_beta_TD_fit.pvalues[1], 3)), fontsize=5, transform=ax6.transAxes)

x_ax6 = np.linspace(7, 11.5, num=500)
y_pred_ax8 = N_0_beta_TD_fit.params[0] + x_ax6 * N_0_beta_TD_fit.params[1]
ax6.plot(10**x_ax6, y_pred_ax8, c='k')

#y1 = list(range(len(taxa_to_keep)))
#latex_labels = [latex_dict[x] for x in taxa_to_keep]
ax1.set_xscale('log', base=10)
ax2.set_xscale('log', base=10)
#ax3.set_xscale('log', basex=10)
ax4.set_xscale('log', base=10)
ax5.set_xscale('log', base=10)
ax6.set_xscale('log', base=10)
#ax7.set_xscale('log', basex=10)
ax8.set_xscale('log', base=10)


ax1.set_xlim([int(8e6), int(7e9)])
ax5.set_xlim([int(8e6), int(7e9)])




ax1.tick_params(axis='x', labelsize=5)
ax2.tick_params(axis='x', labelsize=5)
ax3.tick_params(axis='x', labelsize=5)
ax4.tick_params(axis='x', labelsize=5)
ax5.tick_params(axis='x', labelsize=5)
ax6.tick_params(axis='x', labelsize=5)
ax7.tick_params(axis='x', labelsize=5)
ax8.tick_params(axis='x', labelsize=5)

ax1.tick_params(axis='y', labelsize=4)
ax2.tick_params(axis='y', labelsize=4)
ax3.tick_params(axis='y', labelsize=4)
ax4.tick_params(axis='y', labelsize=4)
ax5.tick_params(axis='y', labelsize=5)
ax6.tick_params(axis='y', labelsize=5)
ax7.tick_params(axis='y', labelsize=5)
ax8.tick_params(axis='y', labelsize=5)


ax5.set_xlabel('Change in population size, ' + r'$\Delta N$', fontsize = 6)
ax6.set_xlabel('Initial death rate, ' + r'$d_{0} \cdot N(0)$', fontsize = 6)
ax7.set_xlabel('Shape parameter, ' + r'$k$', fontsize = 6)
ax8.set_xlabel('Total birth events, ' + r'$n_{births}$', fontsize = 6)

ax1.set_ylabel('Ratio of nonsynonymous\nto synonymous mutations, ' + r'$\frac{dN}{dS}$', fontsize = 6)
ax5.set_ylabel("Tajima's D, " + r'$D_{T}$', fontsize = 6)


fig.savefig(lt.get_path() + '/figs/dn_ds_tajimas_d.pdf', format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
