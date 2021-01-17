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
alpha=1-lt.conf
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
    taxon_color = lt.df_colors.loc[lt.df_colors['strain'] == taxon].Color.to_list()[0]

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


ax.set_xscale('log',base=10)
ax.set_yscale('log',base=10)
ax.set_xlabel('Lag time (hrs.)', fontsize = 16)
ax.set_ylabel('Shape paramter, ' r'$k$', fontsize = 16)
fig.savefig(lt.get_path() + '/figs/lag_shape.pdf', format='pdf',bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
