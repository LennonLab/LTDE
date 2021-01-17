from __future__ import division

import glob, math, os, re, sys
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




#  weibull
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




# innocula 100uL
inoccula = 100
df = pd.read_csv(lt.get_path() + '/data/demography/spoIIE_DC_assay.csv', sep = ',')

#df = pd.read_csv(lt.get_path() + '/data/demography/spo0IIE_assay.csv', sep = ',')

df['N_spores'] = df['HT'] * (1000 / inoccula) * (10 ** (df['dilution_S'] )) * 50 #(mL)
df['N_total'] = df['NT'] * (1000 / inoccula) * (10 ** (df['dilution_V'] )) * 50 #(mL)
df['N_viable'] = df['N_total'] - df['N_spores']
df['days'] = df['hours'] / 24
df = df.sort_values('days')


df_wt = df[(df['strain'] == 'wt')]
df_spoiie = df[(df['strain'] == 'SpoIIE')]


model_wt = log_weibull_model(df_wt.days.values, np.log(df_wt.N_total.values/ df_wt.N_total.values[0]))
result_wt = model_wt.fit(method="lbfgs", disp = False,  bounds= [(0.000001,10000000), (0.00001,100), (0.001, 100)])


model_spoiie = log_weibull_model(df_spoiie.days.values, np.log(df_spoiie.N_total.values/ df_spoiie.N_total.values[0]))
result_spoiie = model_spoiie.fit(method="lbfgs", disp = False,  bounds= [(0.000001,10000000), (0.00001,100), (0.001, 100)])

# percent difference in mean time to death
mttd_wt = lt.get_mttd(result_wt.params[0], result_wt.params[1])
mttd_spoiie = lt.get_mttd(result_spoiie.params[0], result_spoiie.params[1])
percent_difference_mttd = (mttd_spoiie - mttd_wt)/ mttd_wt


k_wt = result_wt.params[1]
d0_wt = result_wt.params[0]

k_spoiie = result_spoiie.params[1]
d0_spoiie = result_spoiie.params[0]



t_ext_wt = ((  -1* np.log(500/df_wt.N_total.values[0] )  )   ** (1/k_wt)) / d0_wt
t_ext_wt = t_ext_wt/365

t_ext_spoiie = ((  -1* np.log(500/df_spoiie.N_total.values[0] )  )   ** (1/k_spoiie)) / d0_spoiie
t_ext_spoiie = t_ext_spoiie/365

percent_difference_t_ext = (t_ext_spoiie - t_ext_wt)/ t_ext_wt



sys.stdout.write("spoIIE mutant, MTTD (days) = %.3f; percent difference = %.3f\n" % (mttd_spoiie, percent_difference_mttd))

sys.stdout.write("spoIIE mutant, time to extinction (years) = %.3f; percent difference = %.3f\n" % (t_ext_spoiie, percent_difference_t_ext))




x_weibull = np.linspace(0, max(df['days'].values))
wt_weibull = df_wt.N_total.values[0] * np.exp( -1 * (( x_weibull* result_wt.params[0]) ** result_wt.params[1]) )
spoiie_weibull = df_spoiie.N_total.values[0] * np.exp( -1 * (( x_weibull* result_spoiie.params[0]) ** result_spoiie.params[1]) )


fig = plt.figure()

taxon_color = lt.df_colors.loc[lt.df_colors['strain'] == 'KBS0812'].Color.to_list()[0]
taxon_color = lt.lighten_color(taxon_color, amount=1.2)

plt.scatter(df_spoiie.days, df_spoiie.N_total, facecolors='none', edgecolors=taxon_color, s=80, lw=2, alpha=0.7, label=lt.latex_dict['KBS0812'] + ' ' + r'$\Delta \mathrm{spoIIE}$',zorder=1)
plt.scatter(df_wt.days, df_wt.N_total, facecolors=taxon_color, edgecolors=taxon_color, s=80, lw=2, alpha=0.7, label=lt.latex_dict['KBS0812'] + ' wt',zorder=1)

plt.plot(x_weibull, wt_weibull, '--', c='k', lw=2, label=lt.latex_dict['KBS0812'] + ' Weibull fit',zorder=2)
plt.plot(x_weibull, spoiie_weibull, ':', c='k', lw=2, label=lt.latex_dict['KBS0812'] + ' ' + r'$\Delta \mathrm{spoIIE}$' + ' Weibull fit',zorder=2)

plt.xlabel('Days, ' + r'$t$', fontsize = 16)
plt.ylabel('Population size, ' + '$N(t)$', fontsize = 16)
plt.yscale('log', base=10)
#plt.ylim(0.5*(10**7), 10**9)
plt.legend(loc='upper right', prop={'size': 8})

#fig.savefig(lt.get_path() + '/figs/spoiie_death_curve.pdf', format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
fig.savefig(lt.get_path() + '/figs/spoiie_DC_death_curve.pdf', format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)

plt.close()
