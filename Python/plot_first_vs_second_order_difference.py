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







to_remove_KBS0711 = [10,11,12]
all_taxa = ['KBS0703', 'ATCC13985', 'ATCC43928', 'KBS0701', 'KBS0702',
            'KBS0705', 'KBS0706', 'KBS0707', 'KBS0710', 'KBS0711',
            'KBS0712', 'KBS0713', 'KBS0714', 'KBS0715', 'KBS0721',
            'KBS0722', 'KBS0724', 'KBS0725', 'KBS0801', 'KBS0802',
            'KBS0812']
df_counts = pd.read_csv(lt.get_path() + '/data/demography/longtermdormancy_20190528_nocomments.csv', sep = ',')
df_counts['Abund'] = (df_counts.Colonies.values+1) * (1000 / df_counts.Inoculum.values) * ( 10 ** df_counts.Dilution.values )
df_counts['Dormstart_date'] = pd.to_datetime(df_counts['Dormstart_date'])
df_counts['Firstread_date'] = pd.to_datetime(df_counts['Firstread_date'])
df_counts['Days'] = df_counts['Firstread_date'] - df_counts['Dormstart_date'] + dt.timedelta(days=1)

fig, ax = plt.subplots(figsize=(4,4))
fig.subplots_adjust(hspace=0.35, wspace=0.35)

for taxon in all_taxa:
    #taxon = 'KBS0812'
    taxon_color = lt.df_colors.loc[lt.df_colors['strain'] == taxon].Color.to_list()[0]
    taxon_genus = lt.df_colors.loc[lt.df_colors['strain'] == taxon].Genus.to_list()[0]
    df_counts_taxon = df_counts.loc[(df_counts['Strain'] == taxon)]
    # ignore 10,11,12
    reps = list(set(df_counts_taxon.Rep.to_list()))
    if taxon == 'KBS0711':
        reps = [j for j in reps if j not in to_remove_KBS0711]
    reps.sort()
    plot_dim = len(reps)/2
    if len(reps) == 4:
        plot_dim_x = 2
        plot_dim_y = 2
        axis_text_font = 8
    elif len(reps) == 6:
        plot_dim_x = 2
        plot_dim_y = 3
        axis_text_font = 6
    else:
        print('Number of reps not recognized')

    for rep in reps:
        df_counts_taxon_rep = df_counts_taxon.loc[(df_counts_taxon['Rep'] == rep)]
        df_counts_taxon_rep = df_counts_taxon_rep.sort_values('Days')

        N_t = df_counts_taxon_rep.Abund.values[2:]
        t = df_counts_taxon_rep.Days.dt.days.values[2:]

        #delta_N_t = N_t[1:] - N_t[:-1]

        #delta_t = t[1:] - t[:-1]

        #delta_N_t_div_t = delta_N_t/delta_t
        ##delta_N_t_div_t = delta_N_t_div_t[np.logical_not(np.isnan(delta_N_t_div_t))]


        #delta_deleta_N_t = delta_N_t[1:] - delta_N_t[:-1]

        #delta_deleta_N_t = delta_N_t_div_t[2:] - 2*delta_N_t_div_t[1:-1] + delta_N_t_div_t[:-2]
        #delta_delta_t = delta_t[1:] - delta_t[:-1]


        #delta_deleta_N_t_div_t = delta_deleta_N_t/(delta_t**2)
        #delta_deleta_N_t_div_t = delta_deleta_N_t_div_t[np.logical_not(np.isnan(delta_deleta_N_t_div_t))]


        #min_delta_N_t = min(delta_N_t_div_t)
        #max_delta_deleta_N_t = max(delta_deleta_N_t_div_t)


        delta_N_t = lt.fd_derivative(N_t, t, n=1, m=2)

        delta_delta_N_t = lt.fd_derivative(N_t, t, n=2, m=2)

        print(np.where(delta_N_t == min(delta_N_t)), np.where(delta_delta_N_t == max(delta_delta_N_t)))

        ax.scatter(abs(min(delta_N_t)), max(delta_delta_N_t), c=taxon_color, marker = 'o', s = 70, \
            linewidth = 0.6, alpha = 0.8, zorder=1, edgecolors='none')

        if taxon == 'KBS0714':

            print(abs(min(delta_N_t)), max(delta_delta_N_t))


        #ax.set_ylim([0.2*min(df_counts_taxon_rep.Abund.values), 4*N_0])

ax.set_xscale('log', base=10)
ax.set_yscale('log', base=10)

ax.set_xlabel('Maximum first-order decrease in ' + r'$N(t)$' + '\nbetween two timepoints, ' + r'$\frac{\left | \Delta N \right |}{ \Delta t}$', fontsize = 12)
ax.set_ylabel('Maximum second-order increase in ' + r'$N(t)$' + '\nbetween two timepoints, ' + r'$\frac{ \Delta^{2} N}{ \Delta t^{2}}$', fontsize = 12)


#fig.text(0.5, 0.02, 'Days, ' + r'$t$', ha='center', fontsize=16)
#fig.text(0.02, 0.5, 'Population size, ' + '$N(t)$', va='center', rotation='vertical', fontsize=16)

#fig.savefig(lt.get_path() + '/figs/taxon_weibull_100/'+taxon+'.png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
fig.savefig(lt.get_path() + '/figs/first_vs_second_order_difference.pdf', format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
#plt.savefig('destination_path.eps', format='eps')

plt.close()
