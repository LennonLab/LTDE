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



fig = plt.figure(figsize = (9, 9))

ax_KBS0714 = plt.subplot2grid((3, 3), (0, 0), colspan=1, rowspan=1)
ax_KBS0703 = plt.subplot2grid((3, 3), (0, 1), colspan=1, rowspan=1)
ax_KBS0812 = plt.subplot2grid((3, 3), (0, 2), colspan=1, rowspan=1)
ax_likelihood = plt.subplot2grid((3, 3), (1, 0), colspan=2, rowspan=1)

ax_KBS0714.set_xlim(-30,1030)
ax_KBS0703.set_xlim(-30,1030)
ax_KBS0812.set_xlim(-30,1030)

ax_KBS0714.text(-0.1, 1.07, 'a', fontsize=13, fontweight='bold', ha='center', va='center', transform=ax_KBS0714.transAxes)
ax_KBS0703.text(-0.1, 1.07, 'b', fontsize=13, fontweight='bold', ha='center', va='center', transform=ax_KBS0703.transAxes)
ax_KBS0812.text(-0.1, 1.07, 'c', fontsize=13, fontweight='bold', ha='center', va='center', transform=ax_KBS0812.transAxes)
ax_likelihood.text(-0.1, 1.07, 'd', fontsize=13, fontweight='bold', ha='center', va='center', transform=ax_likelihood.transAxes)

df_counts = pd.read_csv(lt.get_path() + '/data/demography/longtermdormancy_20190528_nocomments.csv', sep = ',')
df_counts['Abund'] = (df_counts.Colonies.values+1) * (1000 / df_counts.Inoculum.values) * ( 10 ** df_counts.Dilution.values )
df_counts['Dormstart_date'] = pd.to_datetime(df_counts['Dormstart_date'])
df_counts['Firstread_date'] = pd.to_datetime(df_counts['Firstread_date'])
df_counts['Days'] = df_counts['Firstread_date'] - df_counts['Dormstart_date'] + dt.timedelta(days=1)
df_stats = pd.read_csv(lt.get_path() + '/data/demography/weibull_results_clean.csv', sep = ',')

df_counts_KBS0714 = df_counts.loc[((df_counts['Strain'] == 'KBS0714') &  (df_counts['Rep'] == 4))]
df_counts_KBS0703 = df_counts.loc[((df_counts['Strain'] == 'KBS0703') &  (df_counts['Rep'] == 4))]
df_counts_KBS0812 = df_counts.loc[((df_counts['Strain'] == 'KBS0812') &  (df_counts['Rep'] == 4))]

df_stats_KBS0714 = df_stats.loc[((df_stats['strain'] == 'KBS0714') &  (df_stats['rep'] == 4))]
df_stats_KBS0703 = df_stats.loc[((df_stats['strain'] == 'KBS0703') &  (df_stats['rep'] == 4))]
df_stats_KBS0812 = df_stats.loc[((df_stats['strain'] == 'KBS0812') &  (df_stats['rep'] == 4))]

KBS0714_N_0 = df_stats_KBS0714.N_0.to_list()[0]
KBS0703_N_0 = df_stats_KBS0703.N_0.to_list()[0]
KBS0812_N_0 = df_stats_KBS0812.N_0.to_list()[0]

KBS0714_scale = df_stats_KBS0714.beta.to_list()[0]
KBS0703_scale = df_stats_KBS0703.beta.to_list()[0]
KBS0812_scale = df_stats_KBS0812.beta.to_list()[0]

KBS0714_shape = df_stats_KBS0714.alpha.to_list()[0]
KBS0703_shape = df_stats_KBS0703.alpha.to_list()[0]
KBS0812_shape = df_stats_KBS0812.alpha.to_list()[0]

KBS0714_time = list(range(0, max(df_counts_KBS0714.Days.dt.days), 1))
KBS0703_time = list(range(0, max(df_counts_KBS0703.Days.dt.days), 1))
KBS0812_time = list(range(0, max(df_counts_KBS0812.Days.dt.days), 1))

KBS0714_exp_pred = [ KBS0714_N_0* math.exp(-1* (t / KBS0714_scale) ) for t in KBS0714_time]
KBS0703_exp_pred = [ KBS0703_N_0* math.exp(-1* (t / KBS0703_scale) ) for t in KBS0703_time]
KBS0812_exp_pred = [ KBS0812_N_0* math.exp(-1* (t / KBS0812_scale) ) for t in KBS0812_time]

KBS0714_weib_pred = [ KBS0714_N_0* (math.exp(-1* (t / KBS0714_scale)** KBS0714_shape ) )  for t in KBS0714_time]
KBS0703_weib_pred = [ KBS0703_N_0* (math.exp(-1* (t / KBS0703_scale)** KBS0703_shape ) )  for t in KBS0703_time]
KBS0812_weib_pred = [ KBS0812_N_0* (math.exp(-1* (t / KBS0812_scale)** KBS0812_shape ) )  for t in KBS0812_time]

KBS0714_color = lt.df_colors.loc[lt.df_colors['strain'] == 'KBS0714'].Color.to_list()[0]
KBS0703_color = lt.df_colors.loc[lt.df_colors['strain'] == 'KBS0703'].Color.to_list()[0]
KBS0812_color = lt.df_colors.loc[lt.df_colors['strain'] == 'KBS0812'].Color.to_list()[0]


ax_KBS0714.scatter(df_counts_KBS0714.Days.dt.days, df_counts_KBS0714.Abund.values, c=KBS0714_color, marker = 'o', s = 70, \
    linewidth = 0.6, alpha = 0.5, zorder=1, edgecolors='none')
ax_KBS0703.scatter(df_counts_KBS0703.Days.dt.days, df_counts_KBS0703.Abund.values, c=KBS0703_color, marker = 'o', s = 70, \
    linewidth = 0.6, alpha = 0.5, zorder=1, edgecolors='none')
ax_KBS0812.scatter(df_counts_KBS0812.Days.dt.days, df_counts_KBS0812.Abund.values, c=KBS0812_color, marker = 'o', s = 70, \
    linewidth = 0.6, alpha = 0.5, zorder=1, edgecolors='none')

ax_KBS0714.plot(KBS0714_time, KBS0714_exp_pred, zorder=2, c='darkgrey',ls='--', lw=2)
ax_KBS0703.plot(KBS0703_time, KBS0703_exp_pred, zorder=2, c='darkgrey',ls='--', lw=2)
ax_KBS0812.plot(KBS0812_time, KBS0812_exp_pred, zorder=2, c='darkgrey',ls='--', lw=2)

ax_KBS0714.plot(KBS0714_time, KBS0714_weib_pred, zorder=2, c='k',ls='--', lw=2)
ax_KBS0703.plot(KBS0703_time, KBS0703_weib_pred, zorder=2, c='k',ls='--', lw=2)
ax_KBS0812.plot(KBS0812_time, KBS0812_weib_pred, zorder=2, c='k',ls='--', lw=2)

ax_KBS0714.set_ylim([0.2*min(df_counts_KBS0714.Abund.values), 4*KBS0714_N_0])
ax_KBS0703.set_ylim([0.2*min(df_counts_KBS0703.Abund.values), 4*KBS0703_N_0])
ax_KBS0812.set_ylim([0.2*min(df_counts_KBS0812.Abund.values), 4*KBS0812_N_0])

ax_KBS0714.set_yscale('log',base=10)
ax_KBS0703.set_yscale('log',base=10)
ax_KBS0812.set_yscale('log',base=10)

ax_KBS0714.set_title(lt.latex_dict['KBS0714'], fontsize=11)
ax_KBS0703.set_title(lt.latex_dict['KBS0703'], fontsize=11)
ax_KBS0812.set_title(lt.latex_dict['KBS0812'], fontsize=11)

ax_KBS0714.set_xlabel('Days, ' + r'$t$', fontsize=12)
ax_KBS0703.set_xlabel('Days, ' + r'$t$', fontsize=12)
ax_KBS0812.set_xlabel('Days, ' + r'$t$', fontsize=12)

ax_KBS0714.set_ylabel('Population size, ' + '$N_{v}(t)$', fontsize=12)
ax_KBS0703.set_ylabel('Population size, ' + '$N_{v}(t)$', fontsize=12)
ax_KBS0812.set_ylabel('Population size, ' + '$N_{v}(t)$', fontsize=12)

#ax_KBS0714.text(0.65, 0.8,  r'$k=$' + str(round(KBS0714_shape, 2)) , fontsize=9, transform=ax_KBS0714.transAxes)
#ax_KBS0703.text(0.65, 0.8,  r'$k=$' + str(round(KBS0703_shape, 2)) , fontsize=9, transform=ax_KBS0703.transAxes)
#ax_KBS0812.text(0.65, 0.8,  r'$k=$' + str(round(KBS0812_shape, 2)) , fontsize=9, transform=ax_KBS0812.transAxes)


legend_elements_KBS0714 = [Line2D([0], [0], ls='--', color='k', lw=1.5, label='Weibull, ' + r'$k \neq 1$'),
                Line2D([0], [0], ls='--', color='grey', lw=1.5, label= 'Exponen., ' + r'$k = 1$')]

#ax_KBS0714.legend(handles=legend_elements_KBS0714, loc='lower right', fontsize=8)


# plot likelihood
# sort by likelihood

df_stats_mean = df_stats.groupby(['strain'], as_index=False).mean()
df_stats_mean = df_stats_mean.sort_values(by=['LR'])

ax_likelihood.axvline(0, ls='--', c='k', zorder=1)

latex_labels = []
for taxon_idx, taxon in enumerate(df_stats_mean['strain'].values):
    latex_labels.append(lt.latex_dict[taxon])
    taxon_df = df_stats.loc[(df_stats['strain'] == taxon)]
    for index, row in taxon_df.iterrows():
        if row['p.value.BH'] < 0.05:
            color_LR = 'blue'
        else:
            color_LR = 'red'
        ax_likelihood.scatter(row['LR'], taxon_idx, marker = 'o', s = 70, alpha = 0.9, c=color_LR, zorder=2)

legend_elements = [Line2D([0], [0], color = 'none', marker='o', label=r'$P_{\mathrm{BH}}<0.05$',
                    markerfacecolor='b', markeredgecolor='none', markersize=11),
                Line2D([0], [0], marker='o', color='none', label=r'$P_{\mathrm{BH}} \, \nless 0.05$',
                    markerfacecolor='r', markersize=11, markeredgecolor='none', markeredgewidth=2)]

ax_likelihood.legend(handles=legend_elements, loc='lower right')

mttf_taxa_latex = [lt.latex_dict[mttf_taxon] for mttf_taxon in df_stats_mean['strain'].values]
ax_likelihood.yaxis.tick_right()
ax_likelihood.set_yticks(list(range(len(mttf_taxa_latex))))
ax_likelihood.set_yticklabels(mttf_taxa_latex, fontsize=7)

ax_likelihood.set_xlabel('Log-likelihood ratio of the Weibull vs the exponential')



fig.subplots_adjust(wspace=0.35, hspace=0.25)
fig.savefig(lt.get_path() + '/figs/fig1.pdf', format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)

plt.close()
