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
taxa_to_plot = ['ATCC13985', 'KBS0702', 'KBS0707', 'KBS0711', 'KBS0715',
                'KBS0721', 'KBS0722', 'KBS0724', 'KBS0801']


df_colors = pd.read_csv(lt.get_path() + '/data/colors.csv', sep = ',')
conf=0.95

latex_dict = {  'ATCC13985': r'$\mathit{Pseudomonas} \, \mathrm{sp.} \, \mathrm{ATCC13985}$',
                'ATCC43928': r'$\mathit{Pseudomonas} \, \mathrm{sp.} \, \mathrm{ATCC43928}$',
                'KBS0701': r'$\mathit{Pedobacter} \, \mathrm{sp.} \, \mathrm{KBS0701}$',
                'KBS0702': r'$\mathit{Arthrobacter} \, \mathrm{sp.} \, \mathrm{KBS0702}$',
                'KBS0703': r'$\mathit{Arthrobacter} \, \mathrm{sp.} \, \mathrm{KBS0703}$',
                'KBS0705': r'$\mathit{Inquilinus} \, \mathrm{sp.} \, \mathrm{KBS0705}$',
                'KBS0706': r'$\mathit{Mycobacterium} \, \mathrm{sp.} \, \mathrm{KBS0706}$',
                'KBS0707': r'$\mathit{Pseudomonas} \, \mathrm{sp.} \, \mathrm{KBS0707}$',
                'KBS0710': r'$\mathit{Pseudomonas} \, \mathrm{sp.} \, \mathrm{KBS0710}$',
                'KBS0711': r'$\mathit{Janthinobacterium} \, \mathrm{sp.} \, \mathrm{KBS0711}$',
                'KBS0712': r'$\mathit{Variovorax} \, \mathrm{sp.} \, \mathrm{KBS0712}$',
                'KBS0713': r'$\mathit{Yersinia} \, \mathrm{sp.} \, \mathrm{KBS0713}$',
                'KBS0714': r'$\mathit{Micrococcus} \, \mathrm{sp.} \, \mathrm{KBS0714}$',
                'KBS0715': r'$\mathit{Curtobacterium} \, \mathrm{sp.} \, \mathrm{KBS0715}$',
                'KBS0721': r'$\mathit{Flavobacterium} \, \mathrm{sp.} \, \mathrm{KBS0721}$',
                'KBS0722': r'$\mathit{Oerskovia} \, \mathrm{sp.} \, \mathrm{KBS0722}$',
                'KBS0724': r'$\mathit{Rhodococcus} \, \mathrm{sp.} \, \mathrm{KBS0724}$',
                'KBS0725': r'$\mathit{Bradyrhizobium} \, \mathrm{sp.} \, \mathrm{KBS0725}$',
                'KBS0801': r'$\mathit{Burkholderia} \, \mathrm{sp.} \, \mathrm{KBS0801}$',
                'KBS0802': r'$\mathit{Pseudomonas} \, \mathrm{sp.} \, \mathrm{KBS0802}$',
                'KBS0812': r'$\mathit{Bacillus} \, \mathrm{sp.} \, \mathrm{KBS0812}$'
                }





def plot_multiplicity_survival():
    df_par = pd.read_csv(lt.get_path() + '/data/breseq/total_parallelism.txt', sep = '\t' )
    # taxa with at least five genes with a multiplicity greater than one
    # don't include KBS0707 even though it has a significant G score, sicne there's
    # only one gene with a multiplicity greater than 1
    annotate_x = [7.5, 3.75, 5.2, 9]
    fig = plt.figure()
    fig.subplots_adjust(hspace=0.35, wspace=0.35)
    for i in range(0, len(taxa_to_plot)):
        taxon = taxa_to_plot[i]
        df_path = lt.get_path() + '/data/breseq/mult_survival_curves/' + taxon + '.txt'
        df = pd.read_csv(df_path, sep = '\t', index_col=0)
        new_x = [1.0] + df.Mult.tolist() + [df.Mult.tolist()[-1]]
        new_obs_y =[1.0] + df.Obs_fract.tolist() + [ 0.0001]
        new_null_y = [1.0] + df.Null_fract.tolist() + [ 0.0001]

        ax = fig.add_subplot(3, 3, i+1)
        ax.plot(new_x, new_obs_y, '-', c='royalblue', lw=4, alpha = 0.8, zorder=1)
        ax.plot(new_x, new_null_y, '-', c='dimgrey', lw=4, alpha = 0.8, zorder=0)
        ax.set_xlim([0.25, max(new_x)+1])
        x_annotate_position = (max(new_x)+1 - 0.25) * 0.5

        taxon_par = df_par.loc[df_par['Taxon'] == taxon]
        #print(taxon_par.G_score)

        ax.annotate(r'$\Delta \ell= $'+ str(round(float(taxon_par.G_score), 3)), (x_annotate_position, 0.9), fontsize=6)
        if np.log10(float(taxon_par.p_value_BH)) < -3:
            ax.annotate(r'$\mathrm{p_{BH}} = $'+ str('%.2E' % Decimal(float(taxon_par.p_value_BH))), (x_annotate_position, 0.75), fontsize=6)
        else:
            ax.annotate(r'$\mathrm{p_{BH}} = $'+ str(round(float(taxon_par.p_value_BH),3)), (x_annotate_position, 0.75), fontsize=6)

        ax.title.set_text(latex_dict[taxon])
        ax.title.set_fontsize(5.5)
        ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))

        ax.xaxis.set_tick_params(labelsize=6)
        ax.yaxis.set_tick_params(labelsize=6)


    fig.text(0.5, 0.02, 'Gene multiplicity, ' + r'$m$', ha='center', fontsize=16)
    fig.text(0.02, 0.5, 'Fraction mutations ' + r'$\geq m$', va='center', rotation='vertical', fontsize=16)

    fig_name = lt.get_path() + '/figs/mult_survival.pdf'
    fig.savefig(fig_name, format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()


def plot_logpvalue_survival():
    fig = plt.figure()
    fig.subplots_adjust(hspace=0.35, wspace=0.35)
    pstar_dict = pickle.load(open(lt.get_path() + '/data/breseq/p_star.txt', 'rb'))
    for i in range(0, len(taxa_to_plot)):
        taxon = taxa_to_plot[i]
        pstar_i = pstar_dict[taxon][1]
        num_significant_i = pstar_dict[taxon][0] -1
        df_path = lt.get_path() + '/data/breseq/logpvalues/' + taxon + '.txt'
        df = pd.read_csv(df_path, sep = '\t', index_col=0)
        new_x = df.P_value.tolist()
        new_obs_y = df.Obs_num.tolist()
        new_null_y = df.Null_num.tolist()

        ax = fig.add_subplot(3, 3, i+1)

        ax.plot(new_x, new_null_y, '-', c='dimgrey', lw=4, alpha = 0.8, zorder=0)
        ax.plot(new_x, new_obs_y, '-', c='royalblue', lw=4, alpha = 0.8, zorder=1)
        if pstar_i <0:
            y_range = [f[1] for f in list(zip(new_x, new_obs_y)) if f[0] > 0]
            ax.plot([1, 1],[5e-02,max(y_range)],'k-',linewidth=0.5, zorder=2)
            ax.plot([-3,1],[max(y_range), max(y_range)],'k-',linewidth=0.5, zorder=3)
            ax.plot([1], [max(y_range)], c='r', marker='o', zorder=4)
        else:
            ax.plot([pstar_i, pstar_i],[5e-02,num_significant_i],'k-',linewidth=0.5, zorder=2)
            ax.plot([-3,pstar_i],[num_significant_i, num_significant_i],'k-',linewidth=0.5, zorder=3)
            ax.plot([pstar_i], [num_significant_i], c='r', marker='o', zorder=4)

        ax.set_xlim([0.25, max(new_x)+1])

        ax.title.set_text(latex_dict[taxon])
        ax.title.set_fontsize(7.5)

        ax.xaxis.set_tick_params(labelsize=6)
        ax.yaxis.set_tick_params(labelsize=6)
        #pvalue_axis.step(observed_ps/log(10), null_pvalue_survival(observed_ps),'-',label='Expected',color='k')
        #pvalue_axis.step(observed_ps/log(10), observed_pvalue_survival,'b-',label='Observed')

    fig.text(0.5, 0.02, '$-\mathrm{log}_{10}P$', ha='center', fontsize=16)
    fig.text(0.02, 0.5, 'Number of genes', va='center', rotation='vertical', fontsize=16)

    fig_name = lt.get_path() + '/figs/logpvalue_survival.pdf'
    fig.savefig(fig_name, format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()


def distance_decay():
    fig = plt.figure()
    fig.subplots_adjust(hspace=0.35, wspace=0.35)
    for i in range(0, len(taxa_to_plot)):
        taxon = taxa_to_plot[i]
        df_mult_path = lt.get_path() + '/data/breseq/mult_genes/' + taxon + '.txt'
        df_mult = pd.read_csv(df_mult_path, sep = ', ', index_col=0)
        gene_names = df_mult.index.to_list()
        mult_dist = {}
        # get gbff
        # GCF_005937985_2_ASM593798v2_genomic.gbff
        dist_dict = {}
        if taxon == 'KBS0712':
            df_gbff_path = lt.get_path() + '/data/genomes/genomes_ncbi/KBS0712/GCF_006323185_2_ASM632318v2_genomic.gbff'
        elif taxon == 'KBS0711':
            df_gbff_path = lt.get_path() + '/data/genomes/genomes_ncbi/KBS0711/GCF_005937955_2_ASM593795v2_genomic.gbff'
        elif taxon == 'ATCC13985':
            df_gbff_path = lt.get_path() + '/data/genomes/genomes_ncbi/ATCC13985/GCF_006322025_2_ASM632202v2_genomic.gbff'
        elif taxon == 'KBS0715':
            df_gbff_path = lt.get_path() + '/data/genomes/genomes_ncbi/KBS0715/GCF_005938005_2_ASM593800v2_genomic.gbff'
        else:
            continue
        for record in SeqIO.parse(df_gbff_path, "genbank"):
            for feature in record.features:
                if feature.type == 'CDS':
                    if feature.qualifiers['locus_tag'][0] in gene_names:
                        location = str(feature.location)
                        dist_dict[feature.qualifiers['locus_tag'][0]] = int(re.split(r"[:[']+", location)[1])

        dist_df = pd.DataFrame(list(dist_dict.items()), columns=['Name', 'Position'])
        dist_df.index = dist_df.Name
        del dist_df['Name']
        df_merge = df_mult.merge(dist_df, left_index=True, right_index=True)
        ax = fig.add_subplot(2, 2, i+1)
        ax.scatter(df_merge.Position.values, df_merge.Multiplicity.values, c='#175ac6', marker = 'o', s = 70, \
            edgecolors='#244162', linewidth = 0.6, alpha = 0.5, zorder=2)#, edgecolors='none')
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.1e'))

        ax.title.set_text(latex_dict[taxon])

    fig.text(0.5, 0.02, 'Genome position', ha='center', fontsize=16)
    fig.text(0.02, 0.5, 'Gene multiplicity, ' + '$m$', va='center', rotation='vertical', fontsize=16)

    fig.savefig(lt.get_path() + '/figs/dist_mult.png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()


def plot_weib_indiv_taxon():
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
    df_stats = pd.read_csv(lt.get_path() + '/data/demography/weibull_results_clean.csv', sep = ',')
    for taxon in all_taxa:
        print(taxon)
        #taxon = 'KBS0812'
        fig = plt.figure()
        fig.subplots_adjust(hspace=0.35, wspace=0.35)
        taxon_color = df_colors.loc[df_colors['strain'] == taxon].Color.to_list()[0]
        taxon_genus = df_colors.loc[df_colors['strain'] == taxon].Genus.to_list()[0]
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
            df_stats_taxon_rep = df_stats.loc[(df_stats['rep'] == rep) & (df_stats['strain'] == taxon)]
            scale = df_stats_taxon_rep.beta.to_list()[0]
            shape = df_stats_taxon_rep.alpha.to_list()[0]
            N_0 = df_stats_taxon_rep.N_0.to_list()[0]
            last_day = max(df_counts_taxon_rep.Days.dt.days.to_list())
            time = list(range(0, max(df_counts_taxon_rep.Days.dt.days), 1))
            exp_pred = [ N_0* math.exp(-1* (t / scale) ) for t in time]
            weib_pred = [ N_0* (math.exp(-1* (t / scale)** shape ) )  for t in time]
            ax = fig.add_subplot(plot_dim_x, plot_dim_y, rep)
            ax.scatter(df_counts_taxon_rep.Days.dt.days, df_counts_taxon_rep.Abund.values, c=taxon_color, marker = 'o', s = 70, \
                linewidth = 0.6, alpha = 0.5, zorder=1, edgecolors='none')
            ax.plot(time, exp_pred, zorder=2, c='darkgrey',ls='--', lw=2)
            ax.plot(time, weib_pred, zorder=3, c='k',ls='-', lw=2)
            ax.set_ylim([0.2*min(df_counts_taxon_rep.Abund.values), 4*N_0])

            if (scale > 100) or (scale <0.01):
                scale_plot = "{:.2e}".format(scale)
            else:
                scale_plot = round(scale, 2)

            #ax.text(last_day*0.6, 10**(np.log10(4*N_0 - 0.2*min(df_counts_taxon_rep.Abund.values))*0.9), 'Replicate ' + str(rep), fontsize=axis_text_font, transform=ax.transAxes)
            ax.text(0.6, 0.9, 'Replicate ' + str(rep), fontsize=axis_text_font, transform=ax.transAxes)
            ax.text(0.6, 0.77,  r'$d_{0}=$' + str(1/scale_plot), fontsize=axis_text_font, transform=ax.transAxes)
            ax.text(0.6, 0.64,  r'$k=$' + str(round(shape, 2)), fontsize=axis_text_font, transform=ax.transAxes)

            #ax.text(last_day*0.6, 10**(np.log10(4*N_0 - 0.2*min(df_counts_taxon_rep.Abund.values))*0.83),  r'$\lambda=$' + str(scale_plot), fontsize=axis_text_font, transform=ax.transAxes)
            #ax.text(last_day*0.6, 10**(np.log10(4*N_0 - 0.2*min(df_counts_taxon_rep.Abund.values))*0.76),  r'$k=$' + str(round(shape, 2)), fontsize=axis_text_font, transform=ax.transAxes)

            ax.set_yscale('log')

        fig.suptitle(latex_dict[taxon], fontsize=16)

        fig.text(0.5, 0.02, 'Days, ' + r'$t$', ha='center', fontsize=16)
        fig.text(0.02, 0.5, 'Population size, ' + '$N(t)$', va='center', rotation='vertical', fontsize=16)

        #fig.savefig(lt.get_path() + '/figs/taxon_weibull_100/'+taxon+'.png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
        fig.savefig(lt.get_path() + '/figs/taxon_weibull/'+taxon+'.pdf', format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
        #plt.savefig('destination_path.eps', format='eps')

        plt.close()


def mult_syn_nonsyn():
    fig = plt.figure()
    fig.subplots_adjust(hspace=0.35, wspace=0.35)
    for i in range(0, len(taxa_to_plot)):
        taxon = taxa_to_plot[i]
        taxon_color = df_colors.loc[df_colors['strain'] == taxon].Color.to_list()[0]
        df = pd.read_csv(lt.get_path() + '/data/breseq/mult_genes_all/' + taxon + '.txt', sep = ',')

        ax = fig.add_subplot(3, 3, i+1)
        x = df.mult_syn.values
        y = df.mult.values
        x_log10 = np.log10(x)
        y_log10 = np.log10(y)

        ax.scatter(x, y, c=taxon_color, marker = 'o', s = 70, \
            linewidth = 0.6, alpha = 0.5, zorder=1, edgecolors='none')
        min_range = min([min(x), min(y) ] ) * 0.5
        max_range = max([max(x), max(y) ] ) * 2
        ax.set_xlim([min_range, max_range])
        ax.set_ylim([min_range, max_range])
        ax.plot([min_range, max_range], [min_range, max_range], color='darkgrey', linestyle='--', linewidth=2)

        slope, intercept, r_value, p_value, std_err = stats.linregress(x_log10, y_log10)
        ax.set_xscale('log')
        ax.set_yscale('log')

        tt = (slope - 1) / std_err

        p_value_null_1 = 2* (1 - t.cdf(np.abs(tt),df=len(x_log10)-2))

        ax.set_title(latex_dict[taxon])
        ax.title.set_fontsize(7.5)

        ax.xaxis.set_tick_params(labelsize=6)
        ax.yaxis.set_tick_params(labelsize=6)

        x_annotate_position = ((np.log10(ax.get_xlim()[1] - ax.get_xlim()[0])) * 0.85)
        y_annotate_position = ((np.log10(ax.get_ylim()[1] - ax.get_ylim()[0])) * 0.2)
        if p_value_null_1 >= 0.05:
            #ax.annotate(r'$p \nless 0.05$', (x_annotate_position, y_annotate_position), fontsize=6, ha='center', va='center')#, transform = ax.transAxes)
            continue

        x_log10_range = np.linspace( np.log10(min_range), np.log10(max_range), num=1000)
        y_log10_range_pred = np.asarray([ intercept + (x_i*slope) for x_i in  x_log10_range])
        ax.plot(10**x_log10_range, 10**y_log10_range_pred, color='k', linestyle='-', linewidth=2)

        y_log10_pred = np.asarray([intercept + (slope*x_log_10_i) for x_log_10_i in x_log10])
        SSE = sum((y_log10 - y_log10_pred) ** 2)
        N = len(x)
        sd_SSE = np.sqrt( (1/ (N-2)) * SSE)
        sxd=np.sum((x_log10-np.mean(x_log10))**2)

        sx=(x_log10_range-np.mean(x_log10))**2	# x axisr for band
        # Quantile of Student's t distribution for p=1-alpha/2
        alpha=1-conf
        q = stats.t.ppf(1-alpha/2, N-2)
        # Confidence band
        dy = q*sd_SSE*np.sqrt( 1/N + sx/sxd )
        #print(dy)
        # Upper confidence band
        ucb = y_log10_range_pred + dy
        # Lower confidence band
        lcb = y_log10_range_pred - dy

        ax.plot(10**x_log10_range, 10**lcb, color='k', linestyle=':', linewidth=2)
        ax.plot(10**x_log10_range, 10**ucb, color='k', linestyle=':', linewidth=2)

        ax.text(0.65, 0.35, r'$\beta_{1}=$' + str(round(slope, 2 )), fontsize=6, transform=ax.transAxes)
        ax.text(0.65, 0.23, r'$r^{2}=$' + str(round(r_value**2, 2 )), fontsize=6, transform=ax.transAxes)
        ax.text(0.65, 0.11, r'$P<0.05$' , fontsize=6, transform=ax.transAxes)

        # get mean deviation of the data from the linear model given
        # residual sum of squares



    fig.text(0.5, 0.02, 'Synonymous multiplicity', ha='center', fontsize=16)
    fig.text(0.02, 0.5, 'Non-synonymous multiplicity', va='center', rotation='vertical', fontsize=16)
    fig.savefig(lt.get_path() + '/figs/mult_syn.pdf', format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()


def mult_freq():
    fig = plt.figure()
    fig.subplots_adjust(hspace=0.35, wspace=0.35)
    for i in range(0, len(taxa_to_plot)):
        taxon = taxa_to_plot[i]
        taxon_color = df_colors.loc[df_colors['strain'] == taxon].Color.to_list()[0]
        df = pd.read_csv(lt.get_path() + '/data/breseq/mult_genes_all/' + taxon + '.txt', sep = ',')

        ax = fig.add_subplot(3, 3, i+1)
        x = df.mean_freq.values
        y = df.mult.values
        x_log10 = np.log10(x)
        y_log10 = np.log10(y)

        ax.scatter(x, y, c=taxon_color, marker = 'o', s = 70, \
            linewidth = 0.6, alpha = 0.5, zorder=1, edgecolors='none')
        min_range = min([min(x), min(y) ] ) * 0.5
        max_range = max([max(x), max(y) ] ) * 2
        ax.set_xlim([min_range, 1.15])
        ax.set_ylim([min_range, max_range])
        #ax.plot([min_range, max_range], [min_range, max_range], color='darkgrey', linestyle='--', linewidth=2)

        slope, intercept, r_value, p_value, std_err = stats.linregress(x_log10, y_log10)
        ax.set_xscale('log')
        ax.set_yscale('log')

        ax.set_title(latex_dict[taxon])
        ax.title.set_fontsize(7.5)

        ax.xaxis.set_tick_params(labelsize=6)
        ax.yaxis.set_tick_params(labelsize=6)

        print(slope, p_value)

        if p_value >= 0.05:
            continue

        x_log10_range = np.linspace( np.log10(min_range), np.log10(max_range), num=1000)
        y_log10_range_pred = np.asarray([ intercept + (x_i*slope) for x_i in  x_log10_range])
        ax.plot(10**x_log10_range, 10**y_log10_range_pred, color='k', linestyle='-', linewidth=2)

        y_log10_pred = np.asarray([intercept + (slope*x_log_10_i) for x_log_10_i in x_log10])
        SSE = sum((y_log10 - y_log10_pred) ** 2)
        N = len(x)
        sd_SSE = np.sqrt( (1/ (N-2)) * SSE)
        sxd=np.sum((x_log10-np.mean(x_log10))**2)

        sx=(x_log10_range-np.mean(x_log10))**2	# x axisr for band
        # Quantile of Student's t distribution for p=1-alpha/2
        alpha=1-conf
        q = stats.t.ppf(1-alpha/2, N-2)
        # Confidence band
        dy = q*sd_SSE*np.sqrt( 1/N + sx/sxd )
        # Upper confidence band
        ucb = y_log10_range_pred + dy
        # Lower confidence band
        lcb = y_log10_range_pred - dy

        ax.plot(10**x_log10_range, 10**lcb, color='red', linestyle=':', linewidth=2)
        ax.plot(10**x_log10_range, 10**ucb, color='red', linestyle=':', linewidth=2)

        ax.text(0.65, 0.4, r'$\beta_{1}=$' + str(round(slope, 2 )), fontsize=6, transform=ax.transAxes)
        ax.text(0.65, 0.28, r'$r^{2}=$' + str(round(r_value**2, 2 )), fontsize=6, transform=ax.transAxes)
        ax.text(0.65, 0.16, r'$P< 0.05$' , fontsize=6, transform=ax.transAxes)


    fig.text(0.5, 0.02, 'Mean mutation frequency', ha='center', fontsize=16)
    fig.text(0.02, 0.5, 'Multiplicity', va='center', rotation='vertical', fontsize=16)
    fig.savefig(lt.get_path() + '/figs/mult_freq.pdf', format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()



def plot_afs():
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

    for i, taxon in enumerate(taxa_to_plot):
        print(taxon)
        taxon_color = df_colors.loc[df_colors['strain'] == taxon].Color.to_list()[0]
        taxon_reps = [x.split('-')[1] for x in afs_taxa_reps if taxon in x]
        taxon_reps.sort()
        ax = fig.add_subplot(3, 3, i+1)
        taxon_means = []
        taxons_max = []
        for taxon_rep in taxon_reps:
            df = pd.read_csv(lt.get_path() + '/data/breseq/allele_freq_spec/' + taxon + '-'+taxon_rep +'.txt', sep = '\t')
            freqs = df.freq.values
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
            ax.hist(freqs, bins= np.logspace(np.log10(0.1),np.log10(1.0), 50), alpha=0.8,  color = taxon_color, weights=np.zeros_like(freqs) + 1. / len(freqs))
            #print(min(freqs))
            #ax.plot(x_grid_, pdf_, alpha=0.8, lw = 2, color = taxon_color) #, marker='o')
            #ax.plot(10**x_grid_, pdf_, alpha=0.8, lw = 2, color = taxon_color) #, marker='o')
            taxon_means.append(np.mean(freqs))
            taxons_max.append(max(freqs))
            #ax.set_xscale('log')


        ax.title.set_text(latex_dict[taxon])
        ax.title.set_fontsize(7)

        ax.axvline(x=np.mean(taxon_means), color='black', linestyle='--', lw = 1.5, alpha=0.8)
        ax.axvline(x=np.mean(taxons_max), color='red', linestyle='--', lw = 1.5, alpha=0.8)

        ax.tick_params(axis='both', which='major', labelsize=7)
        ax.tick_params(axis='both', which='minor', labelsize=5)

        ax.set_xlim([0.05, 0.82])

        #ax.set_xscale('log', basex=10)


    fig.text(0.5, 0.02, 'Mutation frequency', ha='center', fontsize=16)
    fig.text(0.02, 0.5, 'Frequency', va='center', rotation='vertical', fontsize=16)

    fig.savefig(lt.get_path() + '/figs/afs.pdf', format = 'pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()



def plot_dnds():
    df_taxa = pd.read_csv(lt.get_path() + '/data/breseq/dN_dS_taxa.txt', sep = '\t')
    df_taxa = df_taxa.sort_values(by=['dN_dS_total'])
    taxa_to_keep = df_taxa.Species.to_list()
    df_taxa_samples = pd.read_csv(lt.get_path() + '/data/breseq/genetic_diversity.txt', sep = '\t', index_col=None)
    fig = plt.figure()
    for i, taxon in enumerate(taxa_to_keep):
        #print(taxon)
        #print(df_taxa_samples.loc[df_taxa_samples['Species'] == taxon])
        dn_ds = df_taxa_samples.loc[df_taxa_samples['Species'] == taxon].dn_ds_total.values
        dn_ds_mean = np.mean(dn_ds)
        #print(dn_ds)
        dn_ds_sem = np.std(dn_ds)/ np.sqrt(len(dn_ds))
        taxon_color = df_colors.loc[df_colors['strain'] == taxon].Color.to_list()[0]

        jitter = np.random.normal(i, 0.04, size=len(dn_ds))
        plt.axvline(1, color = 'dimgrey', lw = 2, ls = '--', zorder=1)
        plt.errorbar(dn_ds_mean, i, xerr = dn_ds_sem*2, \
            fmt = 'o', alpha = 0.9, barsabove = True, marker = '.', \
            mfc = 'k', mec = 'k', c = 'k', zorder=3, ms=17)
        plt.scatter(dn_ds, jitter, c=taxon_color, marker = 'o', s = 80, \
            edgecolors='#244162', linewidth = 0, alpha = 1, zorder=2)

        p_value = df_taxa.loc[df_taxa['Species'] == taxon].p_BH.to_list()[0]
        if p_value < 0.05:
            plt.text(dn_ds_mean-0.022, i+0.2, r'$\ast$')


    y1 = list(range(len(taxa_to_keep)))
    latex_labels = [latex_dict[x] for x in taxa_to_keep]

    plt.xlabel('Ratio of nonsynonymous\nto synonymous mutations, ' + r'$\frac{pN}{pS}$', fontsize = 15)

    plt.yticks(y1, latex_labels, rotation=0, fontsize=13)
    fig.savefig(lt.get_path() + '/figs/dn_ds.pdf', format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()



def plot_dn_ds_tajimas_d():
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

        taxon_color = df_colors.loc[df_colors['strain'] == taxon].Color.to_list()[0]

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
    ax1.set_xscale('log', basex=10)
    ax2.set_xscale('log', basex=10)
    #ax3.set_xscale('log', basex=10)
    ax4.set_xscale('log', basex=10)
    ax5.set_xscale('log', basex=10)
    ax6.set_xscale('log', basex=10)
    #ax7.set_xscale('log', basex=10)
    ax8.set_xscale('log', basex=10)


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



def plot_prop_dead_cells():
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
        taxon_color = df_colors.loc[df_colors['strain'] == taxon].Color.to_list()[0]
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
    latex_labels = [latex_dict[x] for x in taxa]


    plt.xlabel('Change in proportion of dead cells', fontsize = 12)
    plt.xlim(-1,1.2)
    plt.yticks(y1, latex_labels, rotation=0)
    fig.savefig(lt.get_path() + '/figs/prop_dead.pdf', format= 'pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()


def plot_irep_shape():
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

        taxon_color = df_colors.loc[df_colors['strain'] == taxon].Color.to_list()[0]
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
    print("D.F. = " + str(len(irep_mean_list)-2))
    print("t = " + str(round((slope-0)/std_err, 4)))
    print("r^2 = " + str(round(r_value**2, 4)))
    print("p-value = " + str(round(p_value, 4)))


def plot_lag_shape():
    df_traits = pd.read_csv(lt.get_path() + '/data/traits/traits.txt', sep = '\t')
    df_traits = df_traits.rename(columns={'Code': 'Species'})
    df_weib = pd.read_csv(lt.get_path() + '/data/demography/weibull_results_clean_species.csv', sep = ',')
    df_merged = df_weib.merge(df_traits, on=['Species'])
    taxa = list(set(df_merged.Species.to_list()))

    #fig = plt.figure()
    fig, ax = plt.subplots(figsize=(4,4))

    ax.axhline(1, color = 'dimgrey', lw = 2, ls = '--', zorder=1)

    x_log10 = np.log10(df_merged['Lag'].values)
    y_log10 = df_merged['alpha.log10'].values
    slope, intercept, r_value, p_value, std_err = stats.linregress(x_log10, y_log10)

    min_range = min([min(10**x_log10), min(10**y_log10) ] ) * 0.5
    max_range = max([max(10**x_log10), max(10**y_log10) ] ) * 2

    x_log10_range = np.linspace( np.log10(min_range), np.log10(max_range), num=1000)
    y_log10_range_pred = np.asarray([ intercept + (x_i*slope) for x_i in  x_log10_range])
    ax.plot(10**x_log10_range, 10**y_log10_range_pred, color='k', linestyle='--', linewidth=2, zorder=2)

    y_log10_pred = np.asarray([intercept + (slope*x_log_10_i) for x_log_10_i in x_log10])
    SSE = sum((y_log10 - y_log10_pred) ** 2)
    N = len(x_log10)
    sd_SSE = np.sqrt( (1/ (N-2)) * SSE)
    sxd=np.sum((x_log10-np.mean(x_log10))**2)

    sx=(x_log10_range-np.mean(x_log10))**2	# x axisr for band
    # Quantile of Student's t distribution for p=1-alpha/2
    alpha=1-conf
    q = stats.t.ppf(1-alpha/2, N-2)
    # Confidence band
    dy = q*sd_SSE*np.sqrt( 1/N + sx/sxd )
    # Upper confidence band
    ucb = y_log10_range_pred + dy
    # Lower confidence band
    lcb = y_log10_range_pred - dy

    ax.plot(10**x_log10_range, 10**lcb, color='red', linestyle=':', linewidth=2, zorder=3)
    ax.plot(10**x_log10_range, 10**ucb, color='red', linestyle=':', linewidth=2, zorder=3)


    for i, taxon in enumerate(taxa):
        taxon_color = df_colors.loc[df_colors['strain'] == taxon].Color.to_list()[0]

        lag = df_merged.loc[df_merged['Species'] == taxon]['Lag'].values
        shape = 10**df_merged.loc[df_merged['Species'] == taxon]['alpha.log10'].values
        ax.scatter(lag, shape, c=taxon_color, marker = 'o', s = 80, \
            edgecolors='#244162', linewidth = 0, alpha = 1, zorder=4)


    print("D.F. = " + str(N-2))
    print("t = " + str(round((slope-0)/std_err, 4)))
    print("r^2 = " + str(round(r_value**2, 4)))
    print("p-value = " + str(round(p_value, 4)))

    ax.set_xlim([0.2,80])
    ax.set_ylim([0.06,1.2])

    ax.text(0.06, 0.8, r'$\beta_{1}=$' + str(round(slope, 3)) , fontsize=9, transform=ax.transAxes)
    ax.text(0.06, 0.74, r'$r^{2}=$' + str(round(r_value**2, 3)) , fontsize=9, transform=ax.transAxes)
    ax.text(0.06, 0.68, r'$P=$' + str(round(p_value, 3)) , fontsize=9, transform=ax.transAxes)


    ax.set_xscale('log',basex=10)
    ax.set_yscale('log',basey=10)
    ax.set_xlabel('Lag time (hrs.)', fontsize = 16)
    ax.set_ylabel('Shape paramter, ' r'$k$', fontsize = 16)
    fig.savefig(lt.get_path() + '/figs/lag_shape.pdf', format='pdf',bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()




def plot_birth_per_death_vs_shape():
    df_gen = pd.read_csv(lt.get_path() + '/data/breseq/genetic_diversity.txt', sep = '\t')
    df_gen = df_gen.rename(columns={'Species': 'strain'})
    df_weib = pd.read_csv(lt.get_path() + '/data/demography/weibull_results_clean.csv', sep = ',')
    df_merged = df_weib.merge(df_gen, on=['strain','rep'])
    taxa = list(set(df_merged.strain.to_list()))

    fig = plt.figure()

    plt.scatter(df_merged.max_N_mut.values, df_merged.alpha.values, color='k', linestyle='--', linewidth=2, zorder=2)

    plt.xlim(min(df_merged.max_N_mut.values), max(df_merged.max_N_mut.values))
    plt.ylim(min(df_merged.alpha.values), 1)

    plt.xscale('log',basex=10)
    plt.yscale('log',basey=10)
    plt.xlabel('births-per-deaths (hrs.)', fontsize = 12)
    plt.ylabel('Shape paramter, ' r'$k$', fontsize = 12)
    fig.savefig(lt.get_path() + '/figs/birth_per_death_vs_shape.png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()


def plot_tajimas_d():
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
        taxon_color = df_colors.loc[df_colors['strain'] == taxon].Color.to_list()[0]

        jitter = np.random.normal(i, 0.04, size=len(x_i))
        plt.axvline(0, color = 'dimgrey', lw = 2, ls = '--', zorder=1)
        plt.errorbar(x_i_mean, i, xerr = x_i_sem*2, \
            fmt = 'o', alpha = 0.9, barsabove = True, marker = '.', \
            mfc = 'k', mec = 'k', c = 'k', zorder=3, ms=17)
        plt.scatter(x_i, jitter, c=taxon_color, marker = 'o', s = 80, \
            edgecolors='#244162', linewidth = 0, alpha = 1, zorder=2)

        p_value = df_taxa.loc[df_taxa['Species'] == taxon].p_BH.to_list()[0]
        if p_value < 0.05:
            plt.text(x_i_mean-0.1, i+0.2, r'$\ast$')


    y1 = list(range(len(taxa_to_keep)))
    latex_labels = [latex_dict[x] for x in taxa_to_keep]

    plt.xlabel("Tajima's D, " + r'$D_{T}$', fontsize = 12)

    plt.yticks(y1, latex_labels, rotation=0)
    fig.savefig(lt.get_path() + '/figs/tajimas_d.pdf', format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()



def plot_bacillus_aa():
    df = pd.read_csv(lt.get_path() + '/data/mzML_Files_forR/Bacillus_AA_Conc_1000d.csv', sep = ',')
    df = df.set_index('Name')
    aa_list = ['Ala', 'Gly', 'Val', 'Leu', 'Ile', 'Pro', 'Met', 'Ser', 'Thr', 'Phe',
            'Asp', 'Glu', 'Orn', 'Lys', 'His', 'Tyr', 'Cys-Cys']
    aa_dict = {'Ala':'Alanine', 'Gly':'Glycine', 'Val':'Valine', 'Leu':'Leucine',
                'Ile': 'Isoleucine', 'Pro':'Proline', 'Met':'Methionine',
                'Ser':'Serine', 'Thr':'Threonine', 'Phe':'Phenylalanine',
                'Asp':'Aspartic Acid', 'Glu':'Glutamic acid', 'Orn':'Arginine',
                'Lys':'Lysine', 'His':'Histidine', 'Tyr':'Tyrosine','Cys-Cys':'Cystine'}
    molar_mass_dict = {'Ala':89.094, 'Gly':75.07, 'Val':117.2, 'Leu':113.2,
                        'Ile':113.2, 'Pro':115.1, 'Met':149.2, 'Ser':105.1,
                        'Thr':119.1, 'Phe':165.2, 'Asp':133.1, 'Glu':147.1,
                        'Orn':174.2, 'Lys':146.2, 'His':155.2, 'Tyr':181.2,
                        'Cys-Cys':240.1}

    bio_reps = ['KBS0812A', 'KBS0812B', 'KBS0812C', 'KBS0812D']
    fig = plt.figure()
    plt.axvline(0, color = 'dimgrey', lw = 2, ls = '--', zorder=1)
    for i, aa in enumerate(aa_list):
        aa_series = df[df.columns[df.columns.to_series().str.contains(aa)]]
        column_name = aa_series.columns.values[0]
        blank = aa_series.loc[ 'std8' , : ].values[0]
        aa_series[column_name] = aa_series[column_name].apply(lambda x: x - blank)
        for bio_rep in bio_reps:
            aa_series_bio_rep = aa_series.loc[ df.index[df.index.to_series().str.contains(bio_rep)] , : ]
            aa_measures = aa_series_bio_rep[column_name].values
            aa_measures = (aa_measures*0.0001) / molar_mass_dict[aa]
            if len(aa_measures) < 3:
                continue
            mean = np.mean(aa_measures)
            sem = np.std(aa_measures) / math.sqrt(len(aa_measures))

            jitter_rep = np.random.normal(i, 0.08)

            plt.errorbar(mean, i, xerr = sem*2, \
                fmt = 'o', alpha = 0.7, barsabove = True, marker = '.', \
                mfc = 'b', mec = 'none', c = 'k', zorder=3, ms=17)

    y1 = list(range(len(aa_list)))
    latex_labels = [aa_dict[x] for x in aa_list]

    plt.xlabel("Blank-corrected molar concentration, mol/L" , fontsize = 12)

    plt.yticks(y1, latex_labels, rotation=0)
    fig.savefig(lt.get_path() + '/figs/bacillus_aa.pdf', format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()







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




def plot_ode():
    ts = np.linspace(0, 1000, 10000)
    N_0 = int(1e9)
    P0 = [N_0, 0, 0]
    colors = ['#FF6347', '#FFA500', '#87CEEB']
    labels = [r'$N(0)\cdot d = 10^{6}$', r'$N(0)\cdot d = 10^{7}$', r'$N(0)\cdot d = 10^{8}$']
    fig = plt.figure()

    for i, d in enumerate([ 0.001, 0.01, 0.1 ]):
        Ps = odeint(dP_dt, P0, ts, args=(d,))
        # Ps = odeint(dP_dt, P0, ts)
        N = Ps[:,0]

        plt.plot(ts, N, "-", c = colors[i], label=labels[i])

    plt.legend(loc='lower left', prop={'size': 8})

    plt.yscale('log',basey=10)

    plt.xlabel('Time, ' + r'$t$' , fontsize = 14)
    plt.ylabel('Number of cells, ' + r'$N(t)$' , fontsize = 14)

    fig.savefig(lt.get_path() + '/figs/ode_example.pdf', format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)




# Modified Gompertz Equation
def log_weibull(t, d_0, k):
    t = np.asarray(t)
    return -1*  ((t * d_0) ** k)


# function to generate confidence intervals based on Fisher Information criteria
def CI_FIC(results):
    # standard errors = square root of the diagnol of a variance-covariance matrix
    ses = np.sqrt(np.absolute(np.diagonal(results.cov_params())))
    cfs = results.params
    lw = cfs - (1.96*ses)
    up = cfs +(1.96*ses)
    return (lw, up)


class log_weibull_model(GenericLikelihoodModel):
    def __init__(self, endog, exog, **kwds):
        super(log_weibull_model, self).__init__(endog, exog, **kwds)
        #print len(exog)

    def nloglikeobs(self, params):
        d_0 = params[0]
        k = params[1]
        z = params[2]
        # probability density function (pdf) is the same as dnorm
        exog_pred = log_weibull(self.endog, d_0 = d_0, k = k)
        # need to flatten the exogenous variable
        LL = -stats.norm.logpdf(self.exog.flatten(), loc=exog_pred, scale=np.exp(z))
        return LL

    def fit(self, start_params=None, maxiter=10000, maxfun=5000, method="bfgs", **kwds):

        if start_params is None:
            d_0_start = 0.01
            k_start = 1
            z_start = 0.8

            start_params = np.array([d_0_start, k_start, z_start])

        return super(log_weibull_model, self).fit(start_params=start_params,
                                maxiter=maxiter, method = method, maxfun=maxfun,
                                **kwds)



def plot_spoIIE():
    # innocula 100uL
    inoccula = 100
    df = pd.read_csv(lt.get_path() + '/data/demography/spo0IIE_assay.csv', sep = ',')
    df['N_spores'] = df['HT'] * (1000 / inoccula) * (10 ** (df['dilution_S'] )) * 50 #(mL)
    df['N_total'] = df['NT'] * (1000 / inoccula) * (10 ** (df['dilution_V'] )) * 50 #(mL)
    df['N_viable'] = df['N_total'] - df['N_spores']
    df['days'] = df['hours'] / 24
    df = df.sort_values('days')


    df_wt = df[(df['strain'] == 'wt')]
    df_spoiie = df[(df['strain'] == 'SpoIIE')]


    model_wt = log_weibull_model(df_wt.days.values, np.log(df_wt.N_total.values/ df_wt.N_total.values[0]))
    result_wt = model_wt.fit(method="lbfgs", disp = False,  bounds= [(0.000001,50), (0.001,100), (0.001, 100)])

    model_spoiie = log_weibull_model(df_spoiie.days.values, np.log(df_spoiie.N_total.values/ df_spoiie.N_total.values[0]))
    result_spoiie = model_spoiie.fit(method="lbfgs", disp = False,  bounds= [(0.000001,50), (0.001,100), (0.001, 100)])

    x_weibull = np.linspace(0, max(df['days'].values))
    wt_weibull = df_wt.N_total.values[0] * np.exp( -1 * (( x_weibull* result_wt.params[0]) ** result_wt.params[1]) )
    spoiie_weibull = df_spoiie.N_total.values[0] * np.exp( -1 * (( x_weibull* result_spoiie.params[0]) ** result_spoiie.params[1]) )


    fig = plt.figure()

    taxon_color = df_colors.loc[df_colors['strain'] == 'KBS0812'].Color.to_list()[0]
    taxon_color = lt.lighten_color(taxon_color, amount=1.2)

    plt.scatter(df_spoiie.days, df_spoiie.N_total, facecolors='none', edgecolors=taxon_color, s=80, lw=2, alpha=0.7, label=latex_dict['KBS0812'] + ' ' + r'$\Delta \mathrm{spoIIE}$',zorder=1)
    plt.scatter(df_wt.days, df_wt.N_total, facecolors=taxon_color, edgecolors=taxon_color, s=80, lw=2, alpha=0.7, label=latex_dict['KBS0812'] + ' wt',zorder=1)

    plt.plot(x_weibull, wt_weibull, '--', c='k', lw=2, label=latex_dict['KBS0812'] + ' Weibull fit',zorder=2)
    plt.plot(x_weibull, spoiie_weibull, ':', c='k', lw=2, label=latex_dict['KBS0812'] + ' ' + r'$\Delta \mathrm{spoIIE}$' + ' Weibull fit',zorder=2)

    plt.xlabel('Days, ' + r'$t$', fontsize = 16)
    plt.ylabel('Population size, ' + '$N(t)$', fontsize = 16)
    plt.yscale('log', basey=10)
    #plt.ylim(0.5*(10**7), 10**9)
    plt.legend(loc='upper right', prop={'size': 8})

    fig.savefig(lt.get_path() + '/figs/spoiie_death_curve.pdf', format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)

    plt.close()



def plot_fig2():

    fig = plt.figure(figsize = (9, 9))

    ax_regression = plt.subplot2grid((4, 4), (0, 0), colspan=2, rowspan=2)
    ax_model = plt.subplot2grid((4, 4), (0, 2), colspan=2, rowspan=2)
    ax_mttd = plt.subplot2grid((4, 4), (2, 0), colspan=3, rowspan=2)

    ax_regression.text(-0.1, 1.07, 'a', fontsize=13, fontweight='bold', ha='center', va='center', transform=ax_regression.transAxes)
    ax_model.text(-0.1, 1.07, 'b', fontsize=13, fontweight='bold', ha='center', va='center', transform=ax_model.transAxes)
    ax_mttd.text(-0.1, 1.07, 'c', fontsize=13, fontweight='bold', ha='center', va='center', transform=ax_mttd.transAxes)

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
        color_taxon = df_colors.loc[df_colors['strain'] == taxon].Color.to_list()[0]
        df_taxon = df_weibull.loc[df_weibull['strain'] == taxon]

        ax_regression.axhline(1, lw=1.5, ls=':',color='grey', zorder=1)
        ax_regression.plot(10**df_CIs.x.values, 10**df_CIs.CI_lower.values, ls=':', c='k',zorder=2)
        ax_regression.plot(10**df_CIs.x.values, 10**df_CIs.CI_upper.values, ls=':', c='k',zorder=2)
        ax_regression.scatter( df_taxon.N_0_beta.values,  df_taxon.alpha.values, c=color_taxon, s=80, alpha=0.8)

    x_log10 = np.log10(np.logspace(6, 14, num=100, endpoint=True, base=10))
    ax_regression.plot(10**x_log10, 10**(model_features_dict['phylom_intercept'] + (x_log10 * model_features_dict['phylom_slope'])), ls='--', c='grey', lw=3)
    ax_regression.plot(10**x_log10, 10**(model_features_dict['lmm_intercept'] + (x_log10 * model_features_dict['lmm_slope'])), ls='--', c='k', lw=3)

    ax_regression.set_xscale('log', basex=10)
    ax_regression.set_yscale('log', basey=10)
    #ax_regression.set_xlim([0.05,1.1])
    ax_regression.set_ylim([0.05,1.2])

    ax_regression.text(0.7, 0.88,  r'$\beta_{1}=$' + str(round(model_features_dict['lmm_slope'], 2)), fontsize=9, transform=ax_regression.transAxes)
    ax_regression.text(0.7, 0.8,  r'$r_{m}^{2}=$' + str(round(model_features_dict['r2_m'], 2)), fontsize=9, transform=ax_regression.transAxes)
    ax_regression.text(0.7, 0.72,  r'$P<10^{-4}$', fontsize=9, transform=ax_regression.transAxes)

    ax_regression.set_xlabel('Initial number of dead cells, ' + r'$d_{0} \cdot N_{v}(0) $', fontsize=12)
    ax_regression.set_ylabel('Degree that growth rate changes, ' + r'$k$', fontsize=13)

    #ax_regression.text(-0.03, 0.5,  "Growth rate decreases time  No chance in growth rate", va='center', rotation='vertical', fontsize=9, transform=ax_regression.transAxes)
    # second figure
    ts = np.linspace(0, 1000, 10000)
    N_0 = int(1e9)
    P0 = [N_0, 0, 0]
    colors = ['#FF6347', '#FFA500', '#87CEEB']
    labels = [r'$d_{0} \cdot N_{v}(0) = 10^{6}$', r'$d_{0} \cdot N_{v}(0) = 10^{7}$', r'$d_{0} \cdot N_{v}(0) = 10^{8}$']

    for i, d in enumerate([ 0.001, 0.01, 0.1 ]):
        Ps = odeint(dP_dt, P0, ts, args=(d,))
        # Ps = odeint(dP_dt, P0, ts)
        N = Ps[:,0]

        ax_model.plot(ts, N, "-", c = colors[i], label=labels[i])

    ax_model.legend(loc='lower left', prop={'size': 9})
    ax_model.set_yscale('log',basey=10)
    ax_model.set_xlabel('Time, ' + r'$t$' , fontsize = 14)
    ax_model.set_ylabel('Number of cells, ' + r'$N_{v}(t)$' , fontsize = 13)

    # mttd figure
    df_weibull_species = pd.read_csv(lt.get_path() + '/data/demography/weibull_results_clean_species.csv', sep=',')
    df_weibull_species = df_weibull_species.sort_values('mttf.log10')
    mttf_taxa = df_weibull_species.Species.values[::-1]
    for taxon_idx, taxon in enumerate(mttf_taxa):
        taxon_color =  df_colors.loc[df_colors['strain'] == taxon].Color.to_list()[0]
        mttf_taxon = df_weibull_species[ df_weibull_species['Species'] ==  taxon].mttf.values[0]
        mttf_log10_taxon = np.log10(mttf_taxon)
        mttf_log10_se_taxon = df_weibull_species[ df_weibull_species['Species'] ==  taxon]['pooled.log10.mttf.se'].values[0]


        l = mlines.Line2D([ 10**(mttf_log10_taxon - (2*mttf_log10_se_taxon)),  10**(mttf_log10_taxon+(2*mttf_log10_se_taxon))], [taxon_idx,taxon_idx],lw=3, c = 'k',zorder=2)
        ax_mttd.add_line(l)
        ax_mttd.scatter(10**mttf_log10_taxon, taxon_idx, marker='o', s = 110, \
                c=taxon_color, alpha=1, zorder=3)

    ax_mttd.axvline( 10**np.mean(df_weibull_species['mttf.log10'].values ), ls = '--', c='grey', lw=2, zorder=1 )

    mttf_taxa_latex = [latex_dict[mttf_taxon] for mttf_taxon in mttf_taxa]
    ax_mttd.yaxis.tick_right()
    ax_mttd.set_yticks(list(range(len(mttf_taxa))))
    ax_mttd.set_yticklabels(mttf_taxa_latex, fontsize=9)
    ax_mttd.set_xlabel('Mean time to death, ' + r'$\bar{T}_{d}$' + ' (days)', fontsize = 18)
    ax_mttd.set_xscale('log', basex=10)

    fig.subplots_adjust(wspace=0.55, hspace=0.45)
    fig.savefig(lt.get_path() + '/figs/fig2.pdf', format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)

    plt.close()



def plot_fig1():

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

    KBS0714_color = df_colors.loc[df_colors['strain'] == 'KBS0714'].Color.to_list()[0]
    KBS0703_color = df_colors.loc[df_colors['strain'] == 'KBS0703'].Color.to_list()[0]
    KBS0812_color = df_colors.loc[df_colors['strain'] == 'KBS0812'].Color.to_list()[0]


    ax_KBS0714.scatter(df_counts_KBS0714.Days.dt.days, df_counts_KBS0714.Abund.values, c=KBS0714_color, marker = 'o', s = 70, \
        linewidth = 0.6, alpha = 0.5, zorder=1, edgecolors='none')
    ax_KBS0703.scatter(df_counts_KBS0703.Days.dt.days, df_counts_KBS0703.Abund.values, c=KBS0703_color, marker = 'o', s = 70, \
        linewidth = 0.6, alpha = 0.5, zorder=1, edgecolors='none')
    ax_KBS0812.scatter(df_counts_KBS0812.Days.dt.days, df_counts_KBS0812.Abund.values, c=KBS0812_color, marker = 'o', s = 70, \
        linewidth = 0.6, alpha = 0.5, zorder=1, edgecolors='none')

    #ax_KBS0714.plot(KBS0714_time, KBS0714_exp_pred, zorder=2, c='darkgrey',ls='--', lw=2)
    #ax_KBS0703.plot(KBS0703_time, KBS0703_exp_pred, zorder=2, c='darkgrey',ls='--', lw=2)
    #ax_KBS0812.plot(KBS0812_time, KBS0812_exp_pred, zorder=2, c='darkgrey',ls='--', lw=2)

    #ax_KBS0714.plot(KBS0714_time, KBS0714_weib_pred, zorder=2, c='k',ls='--', lw=2)
    #ax_KBS0703.plot(KBS0703_time, KBS0703_weib_pred, zorder=2, c='k',ls='--', lw=2)
    #ax_KBS0812.plot(KBS0812_time, KBS0812_weib_pred, zorder=2, c='k',ls='--', lw=2)

    ax_KBS0714.set_ylim([0.2*min(df_counts_KBS0714.Abund.values), 4*KBS0714_N_0])
    ax_KBS0703.set_ylim([0.2*min(df_counts_KBS0703.Abund.values), 4*KBS0703_N_0])
    ax_KBS0812.set_ylim([0.2*min(df_counts_KBS0812.Abund.values), 4*KBS0812_N_0])

    ax_KBS0714.set_yscale('log',basey=10)
    ax_KBS0703.set_yscale('log',basey=10)
    ax_KBS0812.set_yscale('log',basey=10)

    ax_KBS0714.set_title(latex_dict['KBS0714'], fontsize=11)
    ax_KBS0703.set_title(latex_dict['KBS0703'], fontsize=11)
    ax_KBS0812.set_title(latex_dict['KBS0812'], fontsize=11)

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
        latex_labels.append(latex_dict[taxon])
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

    mttf_taxa_latex = [latex_dict[mttf_taxon] for mttf_taxon in df_stats_mean['strain'].values]
    ax_likelihood.yaxis.tick_right()
    ax_likelihood.set_yticks(list(range(len(mttf_taxa_latex))))
    ax_likelihood.set_yticklabels(mttf_taxa_latex, fontsize=7)

    ax_likelihood.set_xlabel('Log-likelihood ratio of the Weibull vs the exponential')



    fig.subplots_adjust(wspace=0.35, hspace=0.25)
    fig.savefig(lt.get_path() + '/figs/fig1.pdf', format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)

    plt.close()




def plot_arthro():


    fig = plt.figure(figsize = (9, 9))

    ax_KBS0703 = plt.subplot2grid((1, 1), (0, 0), colspan=1, rowspan=1)

    ax_KBS0703.set_xlim(-30,1030)

    df_counts = pd.read_csv(lt.get_path() + '/data/demography/longtermdormancy_20190528_nocomments.csv', sep = ',')
    df_counts['Abund'] = (df_counts.Colonies.values+1) * (1000 / df_counts.Inoculum.values) * ( 10 ** df_counts.Dilution.values )
    df_counts['Dormstart_date'] = pd.to_datetime(df_counts['Dormstart_date'])
    df_counts['Firstread_date'] = pd.to_datetime(df_counts['Firstread_date'])
    df_counts['Days'] = df_counts['Firstread_date'] - df_counts['Dormstart_date'] + dt.timedelta(days=1)
    df_stats = pd.read_csv(lt.get_path() + '/data/demography/weibull_results_clean.csv', sep = ',')

    df_counts_KBS0703 = df_counts.loc[((df_counts['Strain'] == 'KBS0703') &  (df_counts['Rep'] == 4))]

    df_stats_KBS0703 = df_stats.loc[((df_stats['strain'] == 'KBS0703') &  (df_stats['rep'] == 4))]

    KBS0703_N_0 = df_stats_KBS0703.N_0.to_list()[0]


    KBS0703_time = list(range(0, max(df_counts_KBS0703.Days.dt.days), 1))

    KBS0714_exp_pred = [ KBS0714_N_0* math.exp(-1* (t / KBS0714_scale) ) for t in KBS0714_time]
    KBS0703_exp_pred = [ KBS0703_N_0* math.exp(-1* (t / KBS0703_scale) ) for t in KBS0703_time]
    KBS0812_exp_pred = [ KBS0812_N_0* math.exp(-1* (t / KBS0812_scale) ) for t in KBS0812_time]

    KBS0714_weib_pred = [ KBS0714_N_0* (math.exp(-1* (t / KBS0714_scale)** KBS0714_shape ) )  for t in KBS0714_time]
    KBS0703_weib_pred = [ KBS0703_N_0* (math.exp(-1* (t / KBS0703_scale)** KBS0703_shape ) )  for t in KBS0703_time]
    KBS0812_weib_pred = [ KBS0812_N_0* (math.exp(-1* (t / KBS0812_scale)** KBS0812_shape ) )  for t in KBS0812_time]

    KBS0714_color = df_colors.loc[df_colors['strain'] == 'KBS0714'].Color.to_list()[0]
    KBS0703_color = df_colors.loc[df_colors['strain'] == 'KBS0703'].Color.to_list()[0]
    KBS0812_color = df_colors.loc[df_colors['strain'] == 'KBS0812'].Color.to_list()[0]


    ax_KBS0714.scatter(df_counts_KBS0714.Days.dt.days, df_counts_KBS0714.Abund.values, c=KBS0714_color, marker = 'o', s = 70, \
        linewidth = 0.6, alpha = 0.5, zorder=1, edgecolors='none')
    ax_KBS0703.scatter(df_counts_KBS0703.Days.dt.days, df_counts_KBS0703.Abund.values, c=KBS0703_color, marker = 'o', s = 70, \
        linewidth = 0.6, alpha = 0.5, zorder=1, edgecolors='none')
    ax_KBS0812.scatter(df_counts_KBS0812.Days.dt.days, df_counts_KBS0812.Abund.values, c=KBS0812_color, marker = 'o', s = 70, \
        linewidth = 0.6, alpha = 0.5, zorder=1, edgecolors='none')

    #ax_KBS0714.plot(KBS0714_time, KBS0714_exp_pred, zorder=2, c='darkgrey',ls='--', lw=2)
    #ax_KBS0703.plot(KBS0703_time, KBS0703_exp_pred, zorder=2, c='darkgrey',ls='--', lw=2)
    #ax_KBS0812.plot(KBS0812_time, KBS0812_exp_pred, zorder=2, c='darkgrey',ls='--', lw=2)

    #ax_KBS0714.plot(KBS0714_time, KBS0714_weib_pred, zorder=2, c='k',ls='--', lw=2)
    #ax_KBS0703.plot(KBS0703_time, KBS0703_weib_pred, zorder=2, c='k',ls='--', lw=2)
    #ax_KBS0812.plot(KBS0812_time, KBS0812_weib_pred, zorder=2, c='k',ls='--', lw=2)

    ax_KBS0714.set_ylim([0.2*min(df_counts_KBS0714.Abund.values), 4*KBS0714_N_0])
    ax_KBS0703.set_ylim([0.2*min(df_counts_KBS0703.Abund.values), 4*KBS0703_N_0])
    ax_KBS0812.set_ylim([0.2*min(df_counts_KBS0812.Abund.values), 4*KBS0812_N_0])

    ax_KBS0714.set_yscale('log',basey=10)
    ax_KBS0703.set_yscale('log',basey=10)
    ax_KBS0812.set_yscale('log',basey=10)

    ax_KBS0714.set_title(latex_dict['KBS0714'], fontsize=11)
    ax_KBS0703.set_title(latex_dict['KBS0703'], fontsize=11)
    ax_KBS0812.set_title(latex_dict['KBS0812'], fontsize=11)

    ax_KBS0714.set_xlabel('Days, ' + r'$t$', fontsize=12)
    ax_KBS0703.set_xlabel('Days, ' + r'$t$', fontsize=12)
    ax_KBS0812.set_xlabel('Days, ' + r'$t$', fontsize=12)

    ax_KBS0714.set_ylabel('Population size, ' + '$N_{v}(t)$', fontsize=12)
    ax_KBS0703.set_ylabel('Population size, ' + '$N_{v}(t)$', fontsize=12)
    ax_KBS0812.set_ylabel('Population size, ' + '$N_{v}(t)$', fontsize=12)





#plot_fig2()

#plot_fig1()
plot_spoIIE()


#plot_dnds()

#plot_tajimas_d()

#plot_bacillus_aa()
#plot_irep_shape()
#plot_logpvalue_survival()
#plot_dnds()

#plot_dn_ds_tajimas_d()

#plot_multiplicity_survival()

#mult_freq()

#plot_afs()
#mult_syn_nonsyn()
#plot_weib_indiv_taxon()

#plot_prop_dead_cells()
#plot_lag_shape()


#plot_weib_indiv_taxon()
