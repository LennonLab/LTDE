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
import matplotlib.gridspec as gridspec

import matplotlib.ticker
import datetime as dt

#from sklearn.model_selection import GridSearchCV
#from sklearn.neighbors import KernelDensity
import statsmodels.stats.multitest as mt
import statsmodels.formula.api as smf

from Bio import SeqIO

from statsmodels.base.model import GenericLikelihoodModel

# in microliters
sample_size = 10
total_volume = 50000

df_colors = lt.df_colors


df_counts = pd.read_csv(lt.get_path() + '/data/demography/longtermdormancy_20190528_nocomments.csv', sep = ',')
df_counts['Abund'] = (df_counts.Colonies.values+1) * (1000 / df_counts.Inoculum.values) * ( 10 ** df_counts.Dilution.values )
df_counts['Dormstart_date'] = pd.to_datetime(df_counts['Dormstart_date'])
df_counts['Firstread_date'] = pd.to_datetime(df_counts['Firstread_date'])
df_counts['Days'] = df_counts['Firstread_date'] - df_counts['Dormstart_date'] + dt.timedelta(days=1)
df_stats = pd.read_csv(lt.get_path() + '/data/demography/weibull_results_clean.csv', sep = ',')

to_remove_KBS0711 = [10,11,12]

df_stats_species = pd.read_csv(lt.get_path() + '/data/demography/weibull_results_clean_species.csv', sep = ',')
df_stats_species = df_stats_species.sort_values(by = 'alpha',  ascending=False)

all_taxa_sorted = df_stats_species.Species.to_list()


all_taxa_nested = [['KBS0703', 'ATCC13985', 'ATCC43928', 'KBS0701', 'KBS0702',
                    'KBS0705', 'KBS0706'], ['KBS0707', 'KBS0710', 'KBS0711',
                    'KBS0712', 'KBS0713', 'KBS0714', 'KBS0715'], ['KBS0721',
                    'KBS0722', 'KBS0724', 'KBS0725', 'KBS0801', 'KBS0802',
                    'KBS0812']]

all_taxa_nested = [all_taxa_sorted[:7], all_taxa_sorted[7:14], all_taxa_sorted[14:]]

#all_taxa_nested = [['KBS0703', 'ATCC13985', 'ATCC43928', 'KBS0701', 'KBS0702',
#                    'KBS0705', 'KBS0706'], ['KBS0707', 'KBS0710', 'KBS0711',
#                    'KBS0712', 'KBS0713', 'KBS0714', 'KBS0715'], ['KBS0721',
#                    'KBS0722', 'KBS0724', 'KBS0725', 'KBS0801', 'KBS0802',
#                    'KBS0812']]


def make_plot():

    for all_taxa_idx, all_taxa in enumerate(all_taxa_nested):

        n_cols = len(all_taxa)

        if all_taxa_idx == 0:
            n_rows = 6
            fig_dims = (12, 8)
        else:
            n_rows = 6
            fig_dims = (12, 8)

        fig = plt.figure(figsize = fig_dims)
        fig.subplots_adjust(hspace=0.35, wspace=0.35)
        gs = gridspec.GridSpec(nrows=n_rows, ncols=n_cols)

        for taxon_idx, taxon in enumerate(all_taxa):

            taxon_color = df_colors.loc[df_colors['strain'] == taxon].Color.to_list()[0]
            taxon_genus = df_colors.loc[df_colors['strain'] == taxon].Genus.to_list()[0]
            df_counts_taxon = df_counts.loc[(df_counts['Strain'] == taxon)]
            # ignore 10,11,12
            reps = list(set(df_counts_taxon.Rep.to_list()))
            if taxon == 'KBS0711':
                reps = [j for j in reps if j not in to_remove_KBS0711]
            reps.sort()

            for rep_idx, rep in enumerate(reps):

                df_counts_taxon_rep = df_counts_taxon.loc[(df_counts_taxon['Rep'] == rep)]
                df_stats_taxon_rep = df_stats.loc[(df_stats['rep'] == rep) & (df_stats['strain'] == taxon)]
                scale = df_stats_taxon_rep.beta.to_list()[0]
                shape = df_stats_taxon_rep.alpha.to_list()[0]
                N_0 = df_stats_taxon_rep.N_0.to_list()[0]
                last_day = max(df_counts_taxon_rep.Days.dt.days.to_list())
                time = list(range(0, max(df_counts_taxon_rep.Days.dt.days), 1))
                exp_pred = [ N_0* math.exp(-1* (t / scale) ) for t in time]
                weib_pred = [ N_0* (math.exp(-1* (t / scale)** shape ) )  for t in time]

                # rows, columns
                ax = fig.add_subplot(gs[rep_idx, taxon_idx])
                #ax = fig.add_subplot(gs[rep_idx, taxon_idx])
                ax.scatter(df_counts_taxon_rep.Days.dt.days, df_counts_taxon_rep.Abund.values, c=taxon_color, marker = 'o', s = 50, \
                    linewidth = 0.6, alpha = 0.5, zorder=1, edgecolors='none')
                ax.plot(time, exp_pred, zorder=2, c='darkgrey',ls='--', lw=2)
                ax.plot(time, weib_pred, zorder=3, c='k',ls='-', lw=2)
                ax.set_ylim([0.2*min(df_counts_taxon_rep.Abund.values), 4*N_0])

                if (scale > 100) or (scale <0.01):
                    scale_plot = "{:.2e}".format(scale)
                else:
                    scale_plot = round(scale, 2)

                #ax.text(last_day*0.6, 10**(np.log10(4*N_0 - 0.2*min(df_counts_taxon_rep.Abund.values))*0.9), 'Replicate ' + str(rep), fontsize=axis_text_font, transform=ax.transAxes)
                ax.text(0.6, 0.9, 'Rep. ' + str(rep), fontsize=6, transform=ax.transAxes)
                #ax.text(0.6, 0.77,  r'$d_{0}=$' + str(1/scale_plot), fontsize=axis_text_font, transform=ax.transAxes)
                #ax.text(0.6, 0.64,  r'$k=$' + str(round(shape, 2)), fontsize=axis_text_font, transform=ax.transAxes)

                #ax.text(last_day*0.6, 10**(np.log10(4*N_0 - 0.2*min(df_counts_taxon_rep.Abund.values))*0.83),  r'$\lambda=$' + str(scale_plot), fontsize=axis_text_font, transform=ax.transAxes)
                #ax.text(last_day*0.6, 10**(np.log10(4*N_0 - 0.2*min(df_counts_taxon_rep.Abund.values))*0.76),  r'$k=$' + str(round(shape, 2)), fontsize=axis_text_font, transform=ax.transAxes)

                ax.set_yscale('log', base=10)
                ax.xaxis.set_tick_params(labelsize=4)
                ax.yaxis.set_tick_params(labelsize=4)

                # labels
                if rep_idx == 0:
                    ax.set_title(lt.latex_dict[taxon],fontsize=5)

                if ( rep_idx == len(reps)-1) :
                    ax.set_xlabel('Days, ' + r'$t$',  fontsize = 6)


                #if rep_idx == len(reps)-2:
                #    if taxon_idx >= len(species_nested_list[rep_idx + 1]):
                #        ax.set_xlabel('Distance between SNVs, $\ell$',  fontsize = 6)

                if taxon_idx == 0:
                    ax.set_ylabel('Population size, ' + '$N(t)$', fontsize = 6)



        #fig.suptitle(latex_dict[taxon], fontsize=16)

        #fig.text(0.5, 0.02, 'Days, ' + r'$t$', ha='center', fontsize=16)
        #fig.text(0.02, 0.5, 'Population size, ' + '$N(t)$', va='center', rotation='vertical', fontsize=16)

        #fig.savefig(lt.get_path() + '/figs/taxon_weibull_100/'+taxon+'.png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
        fig.savefig('%s/figs/survival_curves_all_taxa_%d.pdf' % (lt.get_path(), all_taxa_idx), format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
        #plt.savefig('destination_path.eps', format='eps')
        plt.close()


        #plot_dim = len(reps)/2
        #if len(reps) == 4:
        #    plot_dim_x = 2
        #    plot_dim_y = 2
        #    axis_text_font = 8
        ##elif len(reps) == 6:
        #    plot_dim_x = 2
        #    plot_dim_y = 3
        #    axis_text_font = 6
        #else:
        #    print('Number of reps not recognized')






# get metadata
fraction_sampled_all = []
rep_dict = {}

all_reps = 0
for all_taxa_idx, all_taxa in enumerate(all_taxa_nested):

    for taxon_idx, taxon in enumerate(all_taxa):

        #rep_dict[taxon] = {}

        df_counts_taxon = df_counts.loc[(df_counts['Strain'] == taxon)]
        # ignore 10,11,12
        reps = list(set(df_counts_taxon.Rep.to_list()))

        all_reps += len(reps)

        fraction_sampled_reps = []

        for rep_idx, rep in enumerate(reps):

            df_counts_taxon_rep = df_counts_taxon.loc[(df_counts_taxon['Rep'] == rep)]
            #N_array = df_stats_taxon_rep.N_0.to_list()

            #last_day = df_counts_taxon_rep.Days.dt.days.to_list()

            fraction_sampled = sum(df_counts_taxon_rep.Inoculum.values) / total_volume

            fraction_sampled_reps.append(fraction_sampled)

        rep_dict[taxon] = np.asarray(fraction_sampled_reps)

        fraction_sampled_all.extend(fraction_sampled_reps)


make_plot()

#print(np.mean(fraction_sampled_all))
