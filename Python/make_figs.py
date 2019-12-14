from __future__ import division
import glob, math, os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import ltde_tools as lt
from scipy import stats
from decimal import Decimal
import _pickle as pickle
from matplotlib.ticker import FormatStrFormatter
import matplotlib.ticker
import datetime as dt

from sklearn.model_selection import GridSearchCV
from sklearn.neighbors import KernelDensity
from scipy.stats import t
import statsmodels.stats.multitest as mt

from Bio import SeqIO

# only plot taxa w/ significant g scores and at least 100 mutations
taxa_to_plot = ['ATCC13985', 'KBS0711', 'KBS0712', 'KBS0715']

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

        ax = fig.add_subplot(2, 2, i+1)
        ax.plot(new_x, new_obs_y, '-', c='royalblue', lw=4, alpha = 0.8, zorder=1)
        ax.plot(new_x, new_null_y, '-', c='dimgrey', lw=4, alpha = 0.8, zorder=0)
        ax.set_xlim([0.25, max(new_x)+1])

        taxon_par = df_par.loc[df_par['Taxon'] == taxon]
        #print(taxon_par.G_score)

        ax.annotate(r'$\Delta \ell= $'+ str(round(float(taxon_par.G_score), 3)), (0.6 *max(new_x), 0.9), fontsize=8)
        if np.log10(float(taxon_par.p_value_BH)) < -3:
            ax.annotate(r'$\mathrm{p_{BH}} = $'+ str('%.2E' % Decimal(float(taxon_par.p_value_BH))), (0.6 *max(new_x), 0.75), fontsize=8)
        else:
            ax.annotate(r'$\mathrm{p_{BH}} = $'+ str(round(float(taxon_par.p_value_BH),3)), (0.6 *max(new_x), 0.75), fontsize=8)

        if taxon == 'ATCC13985':
            ax.title.set_text(r'$\mathit{Pseudomonas} \, \mathrm{sp.} \, \mathrm{ATCC13985}$')
        elif taxon == 'KBS0702':
            ax.title.set_text(r'$\mathit{Arthrobacter} \, \mathrm{sp.} \, \mathrm{KBS0702}$')
        elif taxon == 'KBS0711':
            ax.title.set_text(r'$\mathit{Janthinobacterium} \, \mathrm{sp.} \, \mathrm{KBS0711}$')
        elif taxon == 'KBS0712':
            ax.title.set_text(r'$\mathit{Variovorax} \, \mathrm{sp.} \, \mathrm{KBS0712}$')
        elif taxon == 'KBS0713':
            ax.title.set_text(r'$\mathit{Yersinia} \, \mathrm{sp.} \, \mathrm{KBS0713}$')
        elif taxon == 'KBS0715':
            ax.title.set_text(r'$\mathit{Curtobacterium} \, \mathrm{sp.} \, \mathrm{KBS0715}$')
        elif taxon == 'KBS0721':
            ax.title.set_text(r'$\mathit{Flavobacterium} \, \mathrm{sp.} \, \mathrm{KBS0721}$')
        elif taxon == 'KBS0801':
            ax.title.set_text(r'$\mathit{Burkholderia} \, \mathrm{sp.} \, \mathrm{KBS0801}$')
        elif taxon == 'KBS0710':
            ax.title.set_text(r'$\mathit{Pseudomonas} \, \mathrm{sp.} \, \mathrm{KBS0710}$')
        else:
            print("Taxon not recognized")

        ax.title.set_fontsize(8.5)

        ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))

        #ax.set_xticklabels(list(map(int, list(ax.get_xticklabels()))))

    #fig.tight_layout()

    fig.text(0.5, 0.02, 'Gene multiplicity, ' + '$m$', ha='center', fontsize=16)
    fig.text(0.02, 0.5, 'Fraction mutations ' + '$\geq m$', va='center', rotation='vertical', fontsize=16)

    fig_name = lt.get_path() + '/figs/mult_survival.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
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

        ax = fig.add_subplot(2, 2, i+1)

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

        if taxon == 'ATCC13985':
            ax.title.set_text(r'$\mathit{Pseudomonas} \, \mathrm{sp.} \, \mathrm{ATCC13985}$')
        elif taxon == 'KBS0702':
            ax.title.set_text(r'$\mathit{Arthrobacter} \, \mathrm{sp.} \, \mathrm{KBS0702}$')
        elif taxon == 'KBS0711':
            ax.title.set_text(r'$\mathit{Janthinobacterium} \, \mathrm{sp.} \, \mathrm{KBS0711}$')
        elif taxon == 'KBS0712':
            ax.title.set_text(r'$\mathit{Variovorax} \, \mathrm{sp.} \, \mathrm{KBS0712}$')
        elif taxon == 'KBS0713':
            ax.title.set_text(r'$\mathit{Yersinia} \, \mathrm{sp.} \, \mathrm{KBS0713}$')
        elif taxon == 'KBS0715':
            ax.title.set_text(r'$\mathit{Curtobacterium} \, \mathrm{sp.} \, \mathrm{KBS0715}$')
        elif taxon == 'KBS0721':
            ax.title.set_text(r'$\mathit{Flavobacterium} \, \mathrm{sp.} \, \mathrm{KBS0721}$')
        elif taxon == 'KBS0801':
            ax.title.set_text(r'$\mathit{Burkholderia} \, \mathrm{sp.} \, \mathrm{KBS0801}$')
        elif taxon == 'KBS0710':
            ax.title.set_text(r'$\mathit{Pseudomonas} \, \mathrm{sp.} \, \mathrm{KBS0710}$')
        else:
            print("Taxon not recognized")

        ax.title.set_fontsize(8.5)
        #ax.set_xticklabels(tick_labels.astype(int))
        #pvalue_axis.step(observed_ps/log(10), null_pvalue_survival(observed_ps),'-',label='Expected',color='k')

        #pvalue_axis.step(observed_ps/log(10), observed_pvalue_survival,'b-',label='Observed')

    fig.text(0.5, 0.02, '$-\mathrm{log}_{10}P$', ha='center', fontsize=16)
    fig.text(0.02, 0.5, 'Number of genes', va='center', rotation='vertical', fontsize=16)

    fig_name = lt.get_path() + '/figs/logpvalue_survival.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
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

        if taxon == 'ATCC13985':
            ax.title.set_text(r'$\mathit{Pseudomonas} \, \mathrm{sp.} \, \mathrm{ATCC13985}$')
        elif taxon == 'KBS0702':
            ax.title.set_text(r'$\mathit{Arthrobacter} \, \mathrm{sp.} \, \mathrm{KBS0702}$')
        elif taxon == 'KBS0711':
            ax.title.set_text(r'$\mathit{Janthinobacterium} \, \mathrm{sp.} \, \mathrm{KBS0711}$')
        elif taxon == 'KBS0712':
            ax.title.set_text(r'$\mathit{Variovorax} \, \mathrm{sp.} \, \mathrm{KBS0712}$')
        elif taxon == 'KBS0713':
            ax.title.set_text(r'$\mathit{Yersinia} \, \mathrm{sp.} \, \mathrm{KBS0713}$')
        elif taxon == 'KBS0715':
            ax.title.set_text(r'$\mathit{Curtobacterium} \, \mathrm{sp.} \, \mathrm{KBS0715}$')
        elif taxon == 'KBS0721':
            ax.title.set_text(r'$\mathit{Flavobacterium} \, \mathrm{sp.} \, \mathrm{KBS0721}$')
        elif taxon == 'KBS0801':
            ax.title.set_text(r'$\mathit{Burkholderia} \, \mathrm{sp.} \, \mathrm{KBS0801}$')
        elif taxon == 'KBS0710':
            ax.title.set_text(r'$\mathit{Pseudomonas} \, \mathrm{sp.} \, \mathrm{KBS0710}$')
        else:
            print("Taxon not recognized")

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
    df_colors = pd.read_csv(lt.get_path() + '/data/colors.csv', sep = ',')
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

            ax.text(last_day*0.6, 10**(np.log10(4*N_0 - 0.2*min(df_counts_taxon_rep.Abund.values))*0.9), 'Replicate ' + str(rep), fontsize=axis_text_font)
            ax.text(last_day*0.6, 10**(np.log10(4*N_0 - 0.2*min(df_counts_taxon_rep.Abund.values))*0.83),  r'$\lambda=$' + str(scale_plot), fontsize=axis_text_font)
            ax.text(last_day*0.6, 10**(np.log10(4*N_0 - 0.2*min(df_counts_taxon_rep.Abund.values))*0.76),  r'$k=$' + str(round(shape, 2)), fontsize=axis_text_font)

            ax.set_yscale('log')

        if taxon == 'ATCC13985':
            fig.suptitle(r'$\mathit{Pseudomonas} \, \mathrm{sp.} \, \mathrm{ATCC13985}$', fontsize=16)
        elif taxon == 'ATCC43928':
            fig.suptitle(r'$\mathit{Pseudomonas} \, \mathrm{sp.} \, \mathrm{ATCC43928}$', fontsize=16)
        elif taxon == 'KBS0701':
            fig.suptitle(r'$\mathit{Pedobacter} \, \mathrm{sp.} \, \mathrm{KBS0701}$', fontsize=16)
        elif taxon == 'KBS0702':
            fig.suptitle(r'$\mathit{Arthrobacter} \, \mathrm{sp.} \, \mathrm{KBS0702}$', fontsize=16)
        elif taxon == 'KBS0703':
            fig.suptitle(r'$\mathit{Arthrobacter} \, \mathrm{sp.} \, \mathrm{KBS0703}$', fontsize=16)
        elif taxon == 'KBS0705':
            fig.suptitle(r'$\mathit{Inquilinus} \, \mathrm{sp.} \, \mathrm{KBS0705}$', fontsize=16)
        elif taxon == 'KBS0706':
            fig.suptitle(r'$\mathit{Mycobacterium} \, \mathrm{sp.} \, \mathrm{KBS0706}$', fontsize=16)
        elif taxon == 'KBS0707':
            fig.suptitle(r'$\mathit{Pseudomonas} \, \mathrm{sp.} \, \mathrm{KBS0707}$', fontsize=16)
        elif taxon == 'KBS0710':
            fig.suptitle(r'$\mathit{Pseudomonas} \, \mathrm{sp.} \, \mathrm{KBS0710}$', fontsize=16)
        elif taxon == 'KBS0711':
            fig.suptitle(r'$\mathit{Janthinobacterium} \, \mathrm{sp.} \, \mathrm{KBS0711}$', fontsize=16)
        elif taxon == 'KBS0712':
            fig.suptitle(r'$\mathit{Variovorax} \, \mathrm{sp.} \, \mathrm{KBS0712}$', fontsize=16)
        elif taxon == 'KBS0713':
            fig.suptitle(r'$\mathit{Yersinia} \, \mathrm{sp.} \, \mathrm{KBS0713}$', fontsize=16)
        elif taxon == 'KBS0714':
            fig.suptitle(r'$\mathit{Micrococcus} \, \mathrm{sp.} \, \mathrm{KBS0714}$', fontsize=16)
        elif taxon == 'KBS0715':
            fig.suptitle(r'$\mathit{Curtobacterium} \, \mathrm{sp.} \, \mathrm{KBS0715}$', fontsize=16)
        elif taxon == 'KBS0721':
            fig.suptitle(r'$\mathit{Flavobacterium} \, \mathrm{sp.} \, \mathrm{KBS0721}$', fontsize=16)
        elif taxon == 'KBS0722':
            fig.suptitle(r'$\mathit{Oerskovia} \, \mathrm{sp.} \, \mathrm{KBS0722}$', fontsize=16)
        elif taxon == 'KBS0724':
            fig.suptitle(r'$\mathit{Rhodococcus} \, \mathrm{sp.} \, \mathrm{KBS0724}$', fontsize=16)
        elif taxon == 'KBS0725':
            fig.suptitle(r'$\mathit{Bradyrhizobium} \, \mathrm{sp.} \, \mathrm{KBS0725}$', fontsize=16)
        elif taxon == 'KBS0801':
            fig.suptitle(r'$\mathit{Burkholderia} \, \mathrm{sp.} \, \mathrm{KBS0801}$', fontsize=16)
        elif taxon == 'KBS0802':
            fig.suptitle(r'$\mathit{Pseudomonas} \, \mathrm{sp.} \, \mathrm{KBS0802}$', fontsize=16)
        elif taxon == 'KBS0812':
            fig.suptitle(r'$\mathit{Bacillus} \, \mathrm{sp.} \, \mathrm{KBS0812}$', fontsize=16)
        else:
            print("Taxon not recognized")

        fig.text(0.5, 0.02, 'Days, ' + r'$t$', ha='center', fontsize=16)
        fig.text(0.02, 0.5, 'Population size, ' + '$N(t)$', va='center', rotation='vertical', fontsize=16)

        fig.savefig(lt.get_path() + '/figs/taxon_weibull_100/'+taxon+'.png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
        plt.close()


def mult_syn_nonsyn(conf=0.95):
    df_colors = pd.read_csv(lt.get_path() + '/data/colors.csv', sep = ',')
    fig = plt.figure()
    fig.subplots_adjust(hspace=0.35, wspace=0.35)
    for i in range(0, len(taxa_to_plot)):
        taxon = taxa_to_plot[i]
        taxon_color = df_colors.loc[df_colors['strain'] == taxon].Color.to_list()[0]
        df = pd.read_csv(lt.get_path() + '/data/breseq/mult_genes_all/' + taxon + '.txt', sep = ',')

        ax = fig.add_subplot(2, 2, i+1)
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
        print(q)
        # Confidence band
        dy = q*sd_SSE*np.sqrt( 1/N + sx/sxd )
        #print(dy)
        # Upper confidence band
        ucb = y_log10_range_pred + dy
        # Lower confidence band
        lcb = y_log10_range_pred - dy

        ax.plot(10**x_log10_range, 10**lcb, color='red', linestyle=':', linewidth=2)
        ax.plot(10**x_log10_range, 10**ucb, color='red', linestyle=':', linewidth=2)

        # get mean deviation of the data from the linear model given
        # residual sum of squares

        if taxon == 'ATCC13985':
            ax.set_title(r'$\mathit{Pseudomonas} \, \mathrm{sp.} \, \mathrm{ATCC13985}$', fontsize=10)
        elif taxon == 'KBS0711':
            ax.set_title(r'$\mathit{Janthinobacterium} \, \mathrm{sp.} \, \mathrm{KBS0711}$', fontsize=10)
        elif taxon == 'KBS0712':
            ax.set_title(r'$\mathit{Variovorax} \, \mathrm{sp.} \, \mathrm{KBS0712}$', fontsize=10)
        elif taxon == 'KBS0715':
            ax.set_title(r'$\mathit{Curtobacterium} \, \mathrm{sp.} \, \mathrm{KBS0715}$', fontsize=10)
        else:
            continue

    fig.text(0.5, 0.02, 'Synonymous multiplicity', ha='center', fontsize=16)
    fig.text(0.02, 0.5, 'Non-synonymous multiplicity', va='center', rotation='vertical', fontsize=16)
    fig.savefig(lt.get_path() + '/figs/mult_syn.png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()


def mult_freq(conf=0.95):
    df_colors = pd.read_csv(lt.get_path() + '/data/colors.csv', sep = ',')
    fig = plt.figure()
    fig.subplots_adjust(hspace=0.35, wspace=0.35)
    for i in range(0, len(taxa_to_plot)):
        taxon = taxa_to_plot[i]
        taxon_color = df_colors.loc[df_colors['strain'] == taxon].Color.to_list()[0]
        df = pd.read_csv(lt.get_path() + '/data/breseq/mult_genes_all/' + taxon + '.txt', sep = ',')

        ax = fig.add_subplot(2, 2, i+1)
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

        #slope, intercept, r_value, p_value, std_err = stats.linregress(x_log10, y_log10)
        ax.set_xscale('log')
        ax.set_yscale('log')

        #x_log10_range = np.linspace( np.log10(min_range), np.log10(max_range), num=1000)
        #y_log10_range_pred = np.asarray([ intercept + (x_i*slope) for x_i in  x_log10_range])
        #ax.plot(10**x_log10_range, 10**y_log10_range_pred, color='k', linestyle='-', linewidth=2)

        #y_log10_pred = np.asarray([intercept + (slope*x_log_10_i) for x_log_10_i in x_log10])
        #SSE = sum((y_log10 - y_log10_pred) ** 2)
        #N = len(x)
        #sd_SSE = np.sqrt( (1/ (N-2)) * SSE)
        #sxd=np.sum((x_log10-np.mean(x_log10))**2)

        #sx=(x_log10_range-np.mean(x_log10))**2	# x axisr for band
        # Quantile of Student's t distribution for p=1-alpha/2
        #alpha=1-conf
        #q = stats.t.ppf(1-alpha/2, N-2)
        # Confidence band
        #dy = q*sd_SSE*np.sqrt( 1/N + sx/sxd )
        # Upper confidence band
        #ucb = y_log10_range_pred + dy
        # Lower confidence band
        #lcb = y_log10_range_pred - dy

        #ax.plot(10**x_log10_range, 10**lcb, color='red', linestyle=':', linewidth=2)
        #ax.plot(10**x_log10_range, 10**ucb, color='red', linestyle=':', linewidth=2)

        # get mean deviation of the data from the linear model given
        # residual sum of squares

        if taxon == 'ATCC13985':
            ax.set_title(r'$\mathit{Pseudomonas} \, \mathrm{sp.} \, \mathrm{ATCC13985}$', fontsize=10)
        elif taxon == 'KBS0711':
            ax.set_title(r'$\mathit{Janthinobacterium} \, \mathrm{sp.} \, \mathrm{KBS0711}$', fontsize=10)
        elif taxon == 'KBS0712':
            ax.set_title(r'$\mathit{Variovorax} \, \mathrm{sp.} \, \mathrm{KBS0712}$', fontsize=10)
        elif taxon == 'KBS0715':
            ax.set_title(r'$\mathit{Curtobacterium} \, \mathrm{sp.} \, \mathrm{KBS0715}$', fontsize=10)
        else:
            continue

    fig.text(0.5, 0.02, 'Mean mutation frequency', ha='center', fontsize=16)
    fig.text(0.02, 0.5, 'Multiplicity', va='center', rotation='vertical', fontsize=16)
    fig.savefig(lt.get_path() + '/figs/mult_freq.png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()



def plot_afs():
    df_colors = pd.read_csv(lt.get_path() + '/data/colors.csv', sep = ',')
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
    print(str(len(taxa_to_analyze)) +" taxa to analyze")
    fig = plt.figure()
    #taxa_to_analyze = taxa_to_analyze[:2]
    fig.subplots_adjust(hspace=0.35, wspace=0.35)
    for i, taxon in enumerate(taxa_to_analyze):
        print(taxon)
        taxon_color = df_colors.loc[df_colors['strain'] == taxon].Color.to_list()[0]
        taxon_reps = [x.split('-')[1] for x in afs_taxa_reps if taxon in x]
        taxon_reps.sort()
        ax = fig.add_subplot(3, 2, i+1)
        for taxon_rep in taxon_reps:
            df = pd.read_csv(lt.get_path() + '/data/breseq/allele_freq_spec/' + taxon + '-'+taxon_rep +'.txt', sep = '\t')
            #freqs = np.log10(df.freq.values)
            freqs = df.freq.values
            print(freqs)
            #print(len(freqs))
            grid_ = GridSearchCV(KernelDensity(),
                            {'bandwidth': np.linspace(0.001, 1, 50)},
                            cv=2) # 20-fold cross-validation
            #print(freqs[:, None])
            grid_.fit(freqs[:, None])
            #x_grid_ = np.linspace(-2, 0, 1000)
            x_grid_ = np.linspace(0, 1, 1000)
            #print(grid_.best_params_)
            kde_ = grid_.best_estimator_
            pdf_ = np.exp(kde_.score_samples(x_grid_[:, None]))
            #pdf_ = kde_.score_samples(x_grid_[:, None])
            pdf_ = [x / sum(pdf_) for x in pdf_]

            #print(x_grid_)
            ax.plot(x_grid_, pdf_, alpha=0.8, lw = 2, color = taxon_color) #, marker='o')
            #ax.plot(10**x_grid_, pdf_, alpha=0.8, lw = 2, color = taxon_color) #, marker='o')

            #ax.set_xscale('log')
            #ax1.axvline(x=1, color='darkgrey', linestyle='--', lw = 2.5)
            #ax1.axvline(x=np.mean(alpha), color='#1f77b4', linestyle='--', lw = 2.5)

        if taxon == 'ATCC13985':
            ax.title.set_text(r'$\mathit{Pseudomonas} \, \mathrm{sp.} \, \mathrm{ATCC13985}$')
        elif taxon == 'ATCC43928':
            ax.title.set_text(r'$\mathit{Pseudomonas} \, \mathrm{sp.} \, \mathrm{ATCC43928}$')
        elif taxon == 'KBS0701':
            ax.title.set_text(r'$\mathit{Pedobacter} \, \mathrm{sp.} \, \mathrm{KBS0701}$')
        elif taxon == 'KBS0702':
            ax.title.set_text(r'$\mathit{Arthrobacter} \, \mathrm{sp.} \, \mathrm{KBS0702}$')
        elif taxon == 'KBS0703':
            ax.title.set_text(r'$\mathit{Arthrobacter} \, \mathrm{sp.} \, \mathrm{KBS0703}$')
        elif taxon == 'KBS0705':
            ax.title.set_text(r'$\mathit{Inquilinus} \, \mathrm{sp.} \, \mathrm{KBS0705}$')
        elif taxon == 'KBS0706':
            ax.title.set_text(r'$\mathit{Mycobacterium} \, \mathrm{sp.} \, \mathrm{KBS0706}$')
        elif taxon == 'KBS0707':
            ax.title.set_text(r'$\mathit{Pseudomonas} \, \mathrm{sp.} \, \mathrm{KBS0707}$')
        elif taxon == 'KBS0710':
            ax.title.set_text(r'$\mathit{Pseudomonas} \, \mathrm{sp.} \, \mathrm{KBS0710}$')
        elif taxon == 'KBS0711':
            ax.title.set_text(r'$\mathit{Janthinobacterium} \, \mathrm{sp.} \, \mathrm{KBS0711}$')
        elif taxon == 'KBS0712':
            ax.title.set_text(r'$\mathit{Variovorax} \, \mathrm{sp.} \, \mathrm{KBS0712}$')
        elif taxon == 'KBS0713':
            ax.title.set_text(r'$\mathit{Yersinia} \, \mathrm{sp.} \, \mathrm{KBS0713}$')
        elif taxon == 'KBS0714':
            ax.title.set_text(r'$\mathit{Micrococcus} \, \mathrm{sp.} \, \mathrm{KBS0714}$')
        elif taxon == 'KBS0715':
            ax.title.set_text(r'$\mathit{Curtobacterium} \, \mathrm{sp.} \, \mathrm{KBS0715}$')
        elif taxon == 'KBS0721':
            ax.title.set_text(r'$\mathit{Flavobacterium} \, \mathrm{sp.} \, \mathrm{KBS0721}$')
        elif taxon == 'KBS0722':
            ax.title.set_text(r'$\mathit{Oerskovia} \, \mathrm{sp.} \, \mathrm{KBS0722}$')
        elif taxon == 'KBS0724':
            ax.title.set_text(r'$\mathit{Rhodococcus} \, \mathrm{sp.} \, \mathrm{KBS0724}$')
        elif taxon == 'KBS0725':
            ax.title.set_text(r'$\mathit{Bradyrhizobium} \, \mathrm{sp.} \, \mathrm{KBS0725}$')
        elif taxon == 'KBS0801':
            ax.title.set_text(r'$\mathit{Burkholderia} \, \mathrm{sp.} \, \mathrm{KBS0801}$')
        elif taxon == 'KBS0802':
            ax.title.set_text(r'$\mathit{Pseudomonas} \, \mathrm{sp.} \, \mathrm{KBS0802}$')
        elif taxon == 'KBS0812':
            ax.title.set_text(r'$\mathit{Bacillus} \, \mathrm{sp.} \, \mathrm{KBS0812}$')
        else:
            print("Taxon not recognized")

        ax.title.set_fontsize(7)
        ax.tick_params(axis='both', which='major', labelsize=7)
        ax.tick_params(axis='both', which='minor', labelsize=5)


    fig.text(0.5, 0.02, 'Mutation frequency', ha='center', fontsize=16)
    fig.text(0.02, 0.5, 'Frequency', va='center', rotation='vertical', fontsize=16)

    fig.savefig(lt.get_path() + '/figs/afs.png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()



def plot_dnds():
    df_colors = pd.read_csv(lt.get_path() + '/data/colors.csv', sep = ',')
    df_taxa = pd.read_csv(lt.get_path() + '/data/breseq/dN_dS_taxa.txt', sep = '\t')
    df_taxa = df_taxa.sort_values(by=['dN_dS_total'])
    taxa_to_keep = df_taxa.Species.to_list()
    df_taxa_samples = pd.read_csv(lt.get_path() + '/data/breseq/genetic_diversity.txt', sep = '\t', index_col=None)
    fig = plt.figure()
    print(df_taxa_samples)
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
    squad = [r'$\mathit{Pseudomonas} \, \mathrm{sp.} \, \mathrm{ATCC13985}$',
            r'$\mathit{Rhodococcus} \, \mathrm{sp.} \, \mathrm{KBS0724}$',
            r'$\mathit{Variovorax} \, \mathrm{sp.} \, \mathrm{KBS0712}$',
            r'$\mathit{Flavobacterium} \, \mathrm{sp.} \, \mathrm{KBS0721}$',
            r'$\mathit{Burkholderia} \, \mathrm{sp.} \, \mathrm{KBS0801}$',
            r'$\mathit{Curtobacterium} \, \mathrm{sp.} \, \mathrm{KBS0715}$',
            r'$\mathit{Arthrobacter} \, \mathrm{sp.} \, \mathrm{KBS0702}$',
            r'$\mathit{Oerskovia} \, \mathrm{sp.} \, \mathrm{KBS0722}$']

    plt.xlabel('Ratio of nonsynonymous to synonymous mutations, ' + r'$\frac{dN}{dS}$', fontsize = 12)

    plt.yticks(y1, squad, rotation=0)
    fig.savefig(lt.get_path() + '/figs/dn_ds.png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()


def plot_prop_dead_cells():
    df_colors = pd.read_csv(lt.get_path() + '/data/colors.csv', sep = ',')
    df = pd.read_csv(lt.get_path() + '/data/staining.all.new.txt', sep = '\t')
    df = df[df.strain != "KBS0725"]
    taxa = list(set(df.strain.to_list()))
    to_remove = ['KBS0711W', 'KBS0727', 'KBS0714', 'KBS0701']
    taxa = [ x for x in taxa if x not in to_remove]
    df_anc = df.loc[df['hist'] == 'anc']
    df_der = df.loc[df['hist'] == 'der']
    fig = plt.figure()
    plt.axvline(1, color = 'dimgrey', lw = 2, ls = '--', zorder=1)
    # maybe dont include KBS0701
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
        print(taxon)
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


        #tt = (delta_dead_mean)/ delta_dead_sem
        #p_val = t.sf(np.abs(tt), len(delta_dead)-1) # left-tailed one-sided t-test, so use CDF
        #p_values.append(p_val)

    #reject, pvals_corrected, alphacSidak, alphacBonf = mt.multipletests(p_values, alpha=0.05, method='fdr_bh')
    #print(pvals_corrected)
    #for i, taxon in enumerate(taxa_sorted):


    y1 = list(range(len(taxa_sorted)))
    squad = [r'$\mathit{Mycobacterium} \, \mathrm{sp.} \, \mathrm{KBS0706}$',
            r'$\mathit{Arthrobacter} \, \mathrm{sp.} \, \mathrm{KBS0702}$',
            r'$\mathit{Arthrobacter} \, \mathrm{sp.} \, \mathrm{KBS0703}$',
            r'$\mathit{Bacillus} \, \mathrm{sp.} \, \mathrm{KBS0812}$',
            r'$\mathit{Variovorax} \, \mathrm{sp.} \, \mathrm{KBS0712}$',
            r'$\mathit{Pseudomonas} \, \mathrm{sp.} \, \mathrm{KBS0710}$',
            r'$\mathit{Inquilinus} \, \mathrm{sp.} \, \mathrm{KBS0705}$',
            r'$\mathit{Rhodococcus} \, \mathrm{sp.} \, \mathrm{KBS0724}$',
            r'$\mathit{Flavobacterium} \, \mathrm{sp.} \, \mathrm{KBS0721}$',
            r'$\mathit{Yersinia} \, \mathrm{sp.} \, \mathrm{KBS0713}$',
            r'$\mathit{Janthinobacterium} \, \mathrm{sp.} \, \mathrm{KBS0711}$',
            r'$\mathit{Burkholderia} \, \mathrm{sp.} \, \mathrm{KBS0801}$',
            r'$\mathit{Curtobacterium} \, \mathrm{sp.} \, \mathrm{KBS0715}$',
            r'$\mathit{Pseudomonas} \, \mathrm{sp.} \, \mathrm{KBS0802}$',
            r'$\mathit{Pseudomonas} \, \mathrm{sp.} \, \mathrm{KBS0707}$',
            r'$\mathit{Oerskovia} \, \mathrm{sp.} \, \mathrm{KBS0722}$',
            r'$\mathit{Pseudomonas} \, \mathrm{sp.} \, \mathrm{ATCC13985}$',
            r'$\mathit{Pseudomonas} \, \mathrm{sp.} \, \mathrm{ATCC43928}$']


    plt.xlabel('Change in proportion of dead cells', fontsize = 12)
    plt.xlim(-1,1.2)
    plt.yticks(y1, squad, rotation=0)
    fig.savefig(lt.get_path() + '/figs/prop_dead.png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()


#def plot_irep_shape():




#plot_multiplicity_survival()
#plot_logpvalue_survival()
#plot_dnds()
#plot_afs()

#mult_freq()
#plot_afs()
#mult_syn_nonsyn()
#plot_weib_indiv_taxon()

plot_prop_dead_cells()
