import os
import numpy as np
import pandas as pd
import ltde_tools as lt
from statsmodels.base.model import GenericLikelihoodModel
import scipy.stats as stats
import  matplotlib.pyplot as plt

#df_colors = pd.read_csv(lt.get_path() + '/data/colors.csv', sep = ',')

df = pd.read_csv(os.path.expanduser("~/GitHub/LTDE") + '/data/demography/longtermdormancy_20190528_nocomments.csv', sep = ',')
# KBS0721 rep
df['N'] = (df['Colonies']+1) * (1000 / df['Inoculum']) * (10 ** (df['Dilution'] ))

df['Dormstart_date'] =  pd.to_datetime(df['Dormstart_date'], format='%d-%b-%y')
df['Firstread_date'] =  pd.to_datetime(df['Firstread_date'], format='%d-%b-%y')
df['Days'] = df['Firstread_date'].sub(df['Dormstart_date'], axis=0)
df['Days'] = df['Days'].dt.days.astype('int')

taxa = list(set(df.Strain.to_list()))




# Modified Gompertz Equation
def power_law_growth(t, N0, alpha, kappa):
    t = np.asarray(t)
    return N0 / ( (1 - alpha * kappa * t * (N0**alpha) ) ** (1/alpha) )



# function to generate confidence intervals based on Fisher Information criteria
def CI_FIC(results):
    # standard errors = square root of the diagnol of a variance-covariance matrix
    ses = np.sqrt(np.absolute(np.diagonal(results.cov_params())))
    cfs = results.params
    lw = cfs - (1.96*ses)
    up = cfs +(1.96*ses)
    return (lw, up)


class fit_power_law_growth(GenericLikelihoodModel):
    def __init__(self, endog, exog, **kwds):
        super(fit_power_law_growth, self).__init__(endog, exog, **kwds)

    def nloglikeobs(self, params):
        N0 = params[0]
        alpha = params[1]
        kappa = params[2]
        z = params[3]

        # probability density function (pdf) is the same as dnorm
        exog_pred = power_law_growth(self.endog, N0 = N0, alpha = alpha, kappa = kappa)
        # need to flatten the exogenous variable
        LL = -stats.norm.logpdf(self.exog.flatten(), loc=exog_pred, scale=np.exp(z))
        return LL

    def fit(self, start_params=None, maxiter=10000, maxfun=5000, method="bfgs", **kwds):

        if start_params is None:
            N0 = 100000
            alpha_start = 0.1
            kappa_start = 0.0001
            z_start = 0.8

            start_params = np.array([N0, alpha_start, kappa_start, z_start])

        return super(fit_power_law_growth, self).fit(start_params=start_params,
                                maxiter=maxiter, method = method, maxfun=maxfun)





def survival_power_law(t, b, n):
    return -1*b * ( np.asarray(t) **n)


class fit_survival_power_law(GenericLikelihoodModel):
    def __init__(self, endog, exog, **kwds):
        super(fit_power_law_growth, self).__init__(endog, exog, **kwds)

    def nloglikeobs(self, params):
        b = params[0]
        n = params[1]
        z = params[2]

        # probability density function (pdf) is the same as dnorm
        exog_pred = survival_power_law(self.endog, b = b, n = n)
        # need to flatten the exogenous variable
        LL = -stats.norm.logpdf(self.exog.flatten(), loc=exog_pred, scale=np.exp(z))
        return LL

    def fit(self, start_params=None, maxiter=10000, maxfun=5000, method="bfgs", **kwds):

        if start_params is None:
            b = 0.01
            n = 1
            z = 0.1

            start_params = np.array([b, n, z])

        return super(fit_survival_power_law, self).fit(start_params=start_params,
                                maxiter=maxiter, method = method, maxfun=maxfun)







for taxon in taxa:
    if taxon != 'KBS0711':
        continue
    df_taxon = df[(df["Strain"] == taxon)]
    reps = list(set(df_taxon.Rep.to_list()))
    rep = reps[1]

    df_taxon_rep =  df_taxon[(df_taxon["Rep"] == rep)]
    df_taxon_rep = df_taxon_rep.sort_values('Days')


    t = df_taxon_rep.Days.values
    if len(t) < 20:
        continue

    N0 = df_taxon_rep.N.values[0]

    N_t = df_taxon_rep.N.values
    t = df_taxon_rep.Days.values

    S_t = df_taxon_rep.N.values / df_taxon_rep.N.values[0]

    #model = fit_power_law_growth(t, N_t)

    start_params = [N0, 0.1, -0.2, 0.1]

    #result = model.fit(start_params = start_params, method="lbfgs", \
    #                bounds= [(N0,N0), (0.01, 0.9), (-0.1, -0.000000001), (0.001, 10)], \
    #                disp = False)
    model2 = fit_survival_power_law(t, S_t)
    result = model2.fit(start_params = [0.01, 0.1, 0.1], method="bfgs", disp = False)


    print(result.params)

    fig = plt.figure()
    ax = fig.add_subplot(111)


    #N_t_fit = power_law_growth(t, N0, result.params[1], result.params[2])

    #print(power_law_growth(t, N0, result.params[1], result.params[2]))

    #N_t_fit = power_law_growth(t, N0, 2,0.0)

    #print(power_law_growth(t, N0, 0.18,0.001))

    #print(N_t_fit, result.params[1], result.params[2])

    #ax.plot(t, N_t_fit)
    #ax.plot(t, power_law_growth(t, N0, 0.18,-0.001))
    ax.scatter(t, N_t)

    ax.set_yscale('log', basey=10)

    #plt.xlabel(r'$t$', fontsize = 12)
    #plt.ylabel(r'$F(t)$', fontsize = 11)

    fig_name = os.path.expanduser("~/GitHub/LTDE") + '/power_law.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()
