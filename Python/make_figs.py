from __future__ import division
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import ltde_tools as lt
from scipy import stats
from decimal import Decimal


def plot_multiplicity_survival(taxon = 'KBS0711'):
    # taxa with at least five genes with a multiplicity greater than one
    taxa_to_plot = ['ATCC13985', 'KBS0702', 'KBS0711', 'KBS0712', 'KBS0713', 'KBS0715', 'KBS0721', 'KBS0801']


    fig = plt.figure()
    fig.subplots_adjust(hspace=0.45, wspace=0.35)
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
        ax.set_xlim([0.75, 11])

        D, p_value = stats.ks_2samp(df.Null_fract.tolist() , df.Obs_fract.tolist())

        ax.annotate(r'$D= $'+ str(round(D, 3)), (5.2, 0.9), fontsize=8)
        if np.log10(p_value) < -3:
            ax.annotate(r'$p= $'+ str('%.2E' % Decimal(p_value)), (5.2, 0.75), fontsize=8)
        else:
            ax.annotate(r'$p= $'+ str(round(p_value,3)), (5.2, 0.75), fontsize=8)

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
        else:
            print("Taxon not recognized")

        ax.title.set_fontsize(8.5)
    #fig.tight_layout()

    fig.text(0.5, 0.02, 'Gene multiplicity, ' + '$m$', ha='center', fontsize=16)
    fig.text(0.02, 0.5, 'Fraction mutations ' + '$\geq m$', va='center', rotation='vertical', fontsize=16)

    fig_name = lt.get_path() + '/figs/mult_survival.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()





plot_multiplicity_survival()
