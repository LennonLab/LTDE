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






fig = plt.figure()
fig.subplots_adjust(hspace=0.35, wspace=0.35)
for i in range(0, len(lt.taxa_to_plot)):
    taxon = lt.taxa_to_plot[i]
    taxon_color = lt.df_colors.loc[lt.df_colors['strain'] == taxon].Color.to_list()[0]
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

    ax.set_title(lt.latex_dict[taxon])
    ax.title.set_fontsize(7.5)

    ax.xaxis.set_tick_params(labelsize=6)
    ax.yaxis.set_tick_params(labelsize=6)


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
    alpha=1-lt.conf
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
