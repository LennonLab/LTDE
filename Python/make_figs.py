from __future__ import division
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import ltde_tools as lt
from scipy import stats
from decimal import Decimal
import _pickle as pickle
from matplotlib.ticker import FormatStrFormatter

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

        print(taxon_par.G_score)

        ax.annotate(r'$\Delta \ell= $'+ str(round(float(taxon_par.G_score), 3)), (annotate_x[i], 0.9), fontsize=8)
        if np.log10(float(taxon_par.p_value_BH)) < -3:
            ax.annotate(r'$\mathrm{p_{BH}} = $'+ str('%.2E' % Decimal(float(taxon_par.p_value_BH))), (annotate_x[i], 0.75), fontsize=8)
        else:
            ax.annotate(r'$\mathrm{p_{BH}} = $'+ str(round(float(taxon_par.p_value_BH),3)), (annotate_x[i], 0.75), fontsize=8)

        #D, p_value = stats.ks_2samp(df.Null_fract.tolist() , df.Obs_fract.tolist())
        #ax.annotate(r'$D= $'+ str(round(D, 3)), (5.2, 0.9), fontsize=8)
        #if np.log10(p_value) < -3:
        #    ax.annotate(r'$p= $'+ str('%.2E' % Decimal(p_value)), (5.2, 0.75), fontsize=8)
        #else:
        #    ax.annotate(r'$p= $'+ str(round(p_value,3)), (5.2, 0.75), fontsize=8)

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



#plot_multiplicity_survival()
plot_logpvalue_survival()
