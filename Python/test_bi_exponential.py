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


df_colors = pd.read_csv(lt.get_path() + '/data/colors.csv', sep = ',')


import sys
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

#  weibull
def log_weibull(t, d_0, k):
    t = np.asarray(t)
    return np.exp(-1*  ((t * d_0) ** k))


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




def bi_exponential(t, delta, lambda_1, lambda_2):
    t = np.asarray(t)
    return (delta * np.exp(-1*t*lambda_1)) + ((1-delta) * np.exp(-1*t*lambda_2))



class fit_bi_exponential(GenericLikelihoodModel):
    def __init__(self, endog, exog, **kwds):
        super(fit_bi_exponential, self).__init__(endog, exog, **kwds)
        #print len(exog)

    def nloglikeobs(self, params):
        delta = params[0]
        lambda_1 = params[1]
        lambda_2 = params[2]
        z = params[3]
        # probability density function (pdf) is the same as dnorm
        exog_pred = bi_exponential(self.endog, delta = delta, lambda_1 = lambda_1, lambda_2 = lambda_2)
        # need to flatten the exogenous variable
        LL = -stats.norm.logpdf(self.exog.flatten(), loc=exog_pred, scale=np.exp(z))
        return LL

    def fit(self, start_params=None, maxiter=10000, maxfun=5000, method="bfgs", **kwds):

        if start_params is None:
            delta_start = 0.5
            lambda_1_start = 0.03
            lambda_2_start = 0.01
            z_start = 0.8

            start_params = np.array([delta_start, lambda_1_start, lambda_2_start, z_start])

        return super(fit_bi_exponential, self).fit(start_params=start_params,
                                maxiter=maxiter, method = method, maxfun=maxfun,
                                **kwds)





df = pd.read_csv(os.path.expanduser("~/GitHub/LTDE") + '/data/demography/longtermdormancy_20190528_nocomments.csv', sep = ',')
# KBS0721 rep
df['N'] = (df['Colonies']+1) * (1000 / df['Inoculum']) * (10 ** (df['Dilution'] ))

df['Dormstart_date'] =  pd.to_datetime(df['Dormstart_date'], format='%d-%b-%y')
df['Firstread_date'] =  pd.to_datetime(df['Firstread_date'], format='%d-%b-%y')
df['Days'] = df['Firstread_date'].sub(df['Dormstart_date'], axis=0)
df['Days'] = df['Days'].dt.days.astype('int')
#sdf

taxa = list(set(df.Strain.to_list()))

aic_weibull = []
aic_biexponential = []



df_out = open(lt.get_path() + '/data/demography/piecewise_regression.txt', 'w')
df_out.write('\t'.join(['Species', 'rep', 'N0', 'slope1' , 'slope2', 'time_split', 'N_split']) + '\n')

delta_aic_list = []

for taxon in taxa:
    print(taxon)
    if taxon not in df_colors['strain'].values:
        continue
    taxon_color = df_colors.loc[df_colors['strain'] == taxon].Color.to_list()[0]

    df_taxon = df[(df["Strain"] == taxon)]
    reps = list(set(df_taxon.Rep.to_list()))
    for rep in reps:
        df_taxon_rep =  df_taxon[(df_taxon["Rep"] == rep)]
        df_taxon_rep.sort_values('Days')

        t = df_taxon_rep.Days.values
        if len(t) < 20:
            continue
        N = df_taxon_rep.N.values
        N0 = df_taxon_rep.N.values[0]

        proportion = N/N0

        #model_bi_exponential = fit_bi_exponential(t, proportion)
        #result_bi_exponential = model_bi_exponential.fit(method="bfgs", disp = False,  bounds= [(0.99,1), (0.0001,1000), (0.0001,1000), (0.001, 100)])
        log_likelihoods_bi_exponential = []
        for delta in [0.01, 0.1, 1]:
            for lambda_1 in [0.1,1,100]:
                for lambda_2 in [0.1,1,100]:
                    for z in [0.1,1,10]:
                        start_params = [delta, lambda_1, lambda_2, z]

                        #print(t, proportion)
                        model_bi_exponential = fit_bi_exponential(t, proportion)
                        result_bi_exponential = model_bi_exponential.fit(method="lbfgs", disp = False, start_params=start_params, bounds= [(0,1), (0.0001,1000), (0.0001,1000), (0.00001, 100)])
                        ll_bi_exponential = model_bi_exponential.loglike(result_bi_exponential.params)
                        log_likelihoods_bi_exponential.append(ll_bi_exponential)

                        #print(ll_bi_exponential)




        log_likelihoods_weibull = []
        for d_0 in [0.001, 0.01, 0.1, 1, 10, 50, 100, 200]:
            for k in [0.05,0.1,0.5,1]:
                for z in [0.1,1,100]:
                    start_params = [d_0, k, z]

                    model_weibull = log_weibull_model(t, proportion)
                    result_weibull = model_weibull.fit(method="lbfgs", disp = False, start_params=start_params, bounds= [(0.000001,3), (0.00001,100000), (0.0000001, 100)])
                    ll_weibull = model_weibull.loglike(result_weibull.params)
                    #print(ll_weibull)
                    log_likelihoods_weibull.append(ll_weibull)


        ll_bi_exponential = min(log_likelihoods_bi_exponential)
        ll_bi_weibull = min(log_likelihoods_weibull)

        #ll_bi_exponential = model_bi_exponential.loglike(result_bi_exponential.params)

        #n_params_bi_exp = len(result_bi_exponential.params)-1

        #n_params_weibull = len(result_weibull.params)-1

        k_bi_exp = 3
        k_bi_weibull = 2

        aic_bi_exp = 2*k_bi_exp - 2*ll_bi_exponential
        aic_weibull = 2*k_bi_weibull - 2*ll_weibull


        def corrected_aic(aic, n, k):

            return aic + ((2*(k**2) + 2*k) / (n - k - 1))

        aicc_bi_exp = corrected_aic(aic_bi_exp, len(t), k_bi_exp)
        aicc_weibull = corrected_aic(aic_weibull, len(t), k_bi_weibull)


        delta_aic = aicc_weibull - aicc_bi_exp

        delta_aic_list.append(delta_aic)





fig, ax = plt.subplots(figsize=(4,4))

ax.hist(delta_aic_list, bins= 10, color = 'b', density = True)


ax.set_xlabel(r'$cAIC_{\mathrm{Weibull}} - cAIC_{\mathrm{Bi-exponential}}$', fontsize=14)
ax.set_ylabel('Density', fontsize=14)

ax.set_xlim([-4, 2])
ax.axvline(x=0, ls='--', c='k')
#fig.subplots_adjust(wspace=0.3, hspace=0.4)


fig.savefig(lt.get_path() + '/figs/bi_exponential.png', format = 'png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
