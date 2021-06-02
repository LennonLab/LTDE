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






def dP_dt(P, t, d):
    v = 220
    K = 140
    c = 100
    m = 5
    r = 0.005

    uptake=v*(P[2]/(K+P[2]))*P[0] # substrate uptake capacity of viable cells
    cellsMaintained=uptake/(c*m) # how many cells can meet maintenance costs based on substrate uptqke
    cells2die=d*P[0] # how many cells would die if no maintenance costs were met

    dVdt=cellsMaintained-cells2die

    if (cellsMaintained-cells2die) < 0:
        dDdt = cells2die-cellsMaintained
    else:
        dDdt = 0

    dDdt = dDdt - r*P[1]
    dSdt=r*P[1]*c-uptake     #resource

    #dDdt=ifelse((cellsMaintained-cells2die)<0,(cells2die-cellsMaintained),0)-r*D               #cells
    return [dVdt, dDdt, dSdt]






fig = plt.figure(figsize = (9, 17))

ax_regression = plt.subplot2grid((6, 4), (0, 0), colspan=2, rowspan=2)
ax_model = plt.subplot2grid((6, 4), (0, 2), colspan=2, rowspan=2)
ax_mttd = plt.subplot2grid((6, 4), (2, 0), colspan=3, rowspan=2)
ax_ext = plt.subplot2grid((6, 4), (4, 0), colspan=3, rowspan=2)

ax_regression.text(-0.1, 1.07, 'a', fontsize=13, fontweight='bold', ha='center', va='center', transform=ax_regression.transAxes)
ax_model.text(-0.1, 1.07, 'b', fontsize=13, fontweight='bold', ha='center', va='center', transform=ax_model.transAxes)
ax_mttd.text(-0.1, 1.07, 'c', fontsize=13, fontweight='bold', ha='center', va='center', transform=ax_mttd.transAxes)
ax_ext.text(-0.1, 1.07, 'd', fontsize=13, fontweight='bold', ha='center', va='center', transform=ax_ext.transAxes)

df_weibull = pd.read_csv(lt.get_path() + '/data/demography/weibull_results_clean.csv', sep=',')
df_CIs = pd.read_csv(lt.get_path() + '/data/demography/model_CIs.csv', sep=',')

model_features = open(lt.get_path() + '/data/demography/model_features.csv', 'r')
model_features.readline()
model_features_dict = {}
for line in model_features:
    line = line.strip().replace('"', '').split(',')
    model_features_dict[line[0]] = float(line[1])
model_features.close()

taxa = list(set(df_weibull.strain.to_list()))

for taxon in taxa:
    color_taxon = lt.df_colors.loc[lt.df_colors['strain'] == taxon].Color.to_list()[0]
    df_taxon = df_weibull.loc[df_weibull['strain'] == taxon]

    ax_regression.axhline(1, lw=1.5, ls=':',color='grey', zorder=1)
    ax_regression.plot(10**df_CIs.x.values, 10**df_CIs.CI_lower.values, ls=':', c='k',zorder=2)
    ax_regression.plot(10**df_CIs.x.values, 10**df_CIs.CI_upper.values, ls=':', c='k',zorder=2)
    ax_regression.scatter( df_taxon.N_0_beta.values,  df_taxon.alpha.values, c=color_taxon, s=80, alpha=0.8)

x_log10 = np.log10(np.logspace(6, 14, num=100, endpoint=True, base=10))
ax_regression.plot(10**x_log10, 10**(model_features_dict['phylom_intercept'] + (x_log10 * model_features_dict['phylom_slope'])), ls='--', c='grey', lw=3)
ax_regression.plot(10**x_log10, 10**(model_features_dict['lmm_intercept'] + (x_log10 * model_features_dict['lmm_slope'])), ls='--', c='k', lw=3)

ax_regression.set_xscale('log', base=10)
ax_regression.set_yscale('log', base=10)
#ax_regression.set_xlim([0.05,1.1])
ax_regression.set_ylim([0.05,1.2])

ax_regression.text(0.7, 0.88,  r'$\beta_{1}=$' + str(round(model_features_dict['lmm_slope'], 2)), fontsize=11, transform=ax_regression.transAxes)
ax_regression.text(0.7, 0.8,  r'$r_{m}^{2}=$' + str(round(model_features_dict['r2_m'], 2)), fontsize=11, transform=ax_regression.transAxes)
ax_regression.text(0.7, 0.72,  r'$P<10^{-4}$', fontsize=11, transform=ax_regression.transAxes)

#ax_regression.set_xlabel('Initial number of dead cells, ' + r'$d_{0} \cdot N_{v}(0) $', fontsize=12)
ax_regression.set_xlabel('Initial number of dead cells, ' + r'$d_{0} \cdot N(0) $', fontsize=14)

ax_regression.set_ylabel('Degree that growth rate changes, ' + r'$k$', fontsize=14)

#ax_regression.text(-0.03, 0.5,  "Growth rate decreases time  No chance in growth rate", va='center', rotation='vertical', fontsize=9, transform=ax_regression.transAxes)
# second figure
ts = np.linspace(0, 1000, 10000)
N_0 = int(1e9)
P0 = [N_0, 0, 0]
colors = ['#FF6347', '#FFA500', '#87CEEB']
#labels = [r'$d_{0} \cdot N_{v}(0) = 10^{6}$', r'$d_{0} \cdot N_{v}(0) = 10^{7}$', r'$d_{0} \cdot N_{v}(0) = 10^{8}$']
labels = [r'$d_{0} \cdot N(0) = 10^{6}$', r'$d_{0} \cdot N(0) = 10^{7}$', r'$d_{0} \cdot N(0) = 10^{8}$']

for i, d in enumerate([ 0.001, 0.01, 0.1 ]):
    Ps = odeint(dP_dt, P0, ts, args=(d,))
    # Ps = odeint(dP_dt, P0, ts)
    N = Ps[:,0]

    ax_model.plot(ts, N, "-", c = colors[i], label=labels[i])

ax_model.legend(loc='lower left', prop={'size': 9})
ax_model.set_yscale('log',base=10)
ax_model.set_xlabel('Time, ' + r'$t$' , fontsize = 15)
#ax_model.set_ylabel('Number of cells, ' + r'$N_{v}(t)$' , fontsize = 13)
ax_model.set_ylabel('Number of cells, ' + r'$N(t)$' , fontsize = 14)



# mttd figure
df_weibull_species = pd.read_csv(lt.get_path() + '/data/demography/weibull_results_clean_species.csv', sep=',')
df_weibull_species = df_weibull_species.sort_values('mttf.log10')
mttf_taxa = df_weibull_species.Species.values[::-1]
for taxon_idx, taxon in enumerate(mttf_taxa):
    taxon_color =  lt.df_colors.loc[lt.df_colors['strain'] == taxon].Color.to_list()[0]

    #mttf_taxon = df_weibull_species[ df_weibull_species['Species'] ==  taxon].mttf.values[0]
    mttf_taxon = df_weibull_species[ df_weibull_species['Species'] ==  taxon]['mttf.log10'].values[0]

    #mttf_log10_taxon = np.log10(mttf_taxon)
    mttf_log10_taxon = mttf_taxon
    mttf_log10_se_taxon = df_weibull_species[ df_weibull_species['Species'] ==  taxon]['pooled.log10.mttf.se'].values[0]


    l = mlines.Line2D([ 10**(mttf_log10_taxon - (2*mttf_log10_se_taxon)),  10**(mttf_log10_taxon+(2*mttf_log10_se_taxon))], [taxon_idx,taxon_idx],lw=3, c = 'k',zorder=2)
    ax_mttd.add_line(l)
    ax_mttd.scatter(10**mttf_log10_taxon, taxon_idx, marker='o', s = 80, \
            c=taxon_color, alpha=1, zorder=3)

ax_mttd.axvline( 10**np.mean(df_weibull_species['mttf.log10'].values ), ls = '--', c='grey', lw=2, zorder=1 )

mttf_taxa_latex = [lt.latex_dict[mttf_taxon] for mttf_taxon in mttf_taxa]
ax_mttd.yaxis.tick_right()
ax_mttd.set_yticks(list(range(len(mttf_taxa))))
ax_mttd.set_yticklabels(mttf_taxa_latex, fontsize=12)
ax_mttd.set_xlabel('Mean time to death, ' + r'$\bar{T}_{d}$' + ' (days)', fontsize = 18)
ax_mttd.set_xscale('log', base=10)


# time to extinction figure
df_weibull_species = df_weibull_species.sort_values('T_ext.log10')
t_ext_taxa = df_weibull_species.Species.values[::-1]
for taxon_idx, taxon in enumerate(t_ext_taxa):
    taxon_color =  lt.df_colors.loc[lt.df_colors['strain'] == taxon].Color.to_list()[0]
    t_ext_taxon = df_weibull_species[ df_weibull_species['Species'] ==  taxon]['T_ext.log10'].values[0]
    t_ext_log10_taxon = t_ext_taxon

    t_ext_log10_se_taxon = df_weibull_species[ df_weibull_species['Species'] ==  taxon]['pooled.log10.T_ext.se'].values[0]
    l = mlines.Line2D([ (10**(t_ext_log10_taxon - (2*t_ext_log10_se_taxon)))/365,  (10**(t_ext_log10_taxon+(2*t_ext_log10_se_taxon)))/365], [taxon_idx,taxon_idx],lw=3, c = 'k',zorder=2)
    ax_ext.add_line(l)
    ax_ext.scatter((10**t_ext_log10_taxon)/365, taxon_idx, marker='o', s = 80, \
            c=taxon_color, alpha=1, zorder=3)


ax_ext.axvline( (10**np.mean(df_weibull_species['T_ext.log10'].values ))/365, ls = '--', c='grey', lw=2, zorder=1 )

t_ext_taxa_latex = [lt.latex_dict[mttf_taxon] for mttf_taxon in t_ext_taxa]
ax_ext.yaxis.tick_right()
ax_ext.set_yticks(list(range(len(t_ext_taxa))))
ax_ext.set_yticklabels(t_ext_taxa_latex, fontsize=12)
ax_ext.set_xlabel('Time to extinction, ' + r'$T_{ext}$' + ' (years)', fontsize = 18)
ax_ext.set_xscale('log', base=10)


fig.subplots_adjust(wspace=0.55, hspace=0.6)
fig.savefig(lt.get_path() + '/figs/fig2.pdf', format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)

plt.close()
