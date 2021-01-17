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





df_par = pd.read_csv(lt.get_path() + '/data/breseq/total_parallelism.txt', sep = '\t' )
# taxa with at least five genes with a multiplicity greater than one
# don't include KBS0707 even though it has a significant G score, sicne there's
# only one gene with a multiplicity greater than 1
annotate_x = [7.5, 3.75, 5.2, 9]
fig = plt.figure()
fig.subplots_adjust(hspace=0.35, wspace=0.35)
for i in range(0, len(lt.taxa_to_plot)):
    taxon = lt.taxa_to_plot[i]
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

    ax.title.set_text(lt.latex_dict[taxon])
    ax.title.set_fontsize(5.5)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))

    ax.xaxis.set_tick_params(labelsize=6)
    ax.yaxis.set_tick_params(labelsize=6)


fig.text(0.5, 0.02, 'Gene multiplicity, ' + r'$m$', ha='center', fontsize=16)
fig.text(0.02, 0.5, 'Fraction mutations ' + r'$\geq m$', va='center', rotation='vertical', fontsize=16)

fig_name = lt.get_path() + '/figs/mult_survival.pdf'
fig.savefig(fig_name, format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
