from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import statsmodels.api as sm
import ltde_tools as lt

from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV
from scipy import stats
from scipy.optimize import brentq

#df = pd.read_csv(lt.get_path() + '/data/demography/weibull_results_clean.csv', sep = ',')
# remove bacillus
#df = df[df.strain != 'KBS0812']
#df = df[df['beta']>1]
#df_diversity = pd.read_csv(lt.get_path() + '/data/breseq/genetic_diversity.txt', sep = '\t')
#df_merged = df.merge(df_diversity, on=['strain', 'rep'])


def get_random_joint_kde_samples(size):
    print(df_merged.beta.values)
    joint_data = np.log(df_merged[['beta', 'N_0']].values)
    kde = sm.nonparametric.KDEMultivariate(data=joint_data, var_type='cc', bw='normal_reference')
    n, d = kde.data.shape
    indices = np.random.randint(0, n, size)
    cov = np.diag(kde.bw)**2
    means = kde.data[indices, :]
    norm = np.random.multivariate_normal(np.zeros(d), cov, size)
    samples = np.exp(np.transpose(means + norm))
    return samples[0], samples[1]

def sample_divisions(size):
    divisions_np = df_merged.binary_divisions.values
    divisions_np = divisions_np.astype('double')
    kde_birth_events = sm.nonparametric.KDEMultivariate(np.log(divisions_np), var_type='c',  bw='normal_reference')
    # sample
    samples = []
    for x in range(size):
        u = np.random.random()
        # 1-d root-finding
        def func(x):
            return kde_birth_events.cdf([x]) - u
        sample_x = brentq(func, -99999999, 99999999)  # read brentq-docs about these constants
                                                      # constants need to be sign-changing for the function
        samples.append(np.exp(sample_x))
    return samples


def birth_rate(N, S, c, K = 500, alpha=160, V_m=0.3):
    return (1/c) * ((V_m*S) / (K+S)) * (S / (S + alpha*N))


def birth_rate_no_N_depen(S, c, K = 500, V_m=0.3):
    return (1/c) * ((V_m*S) / (K+S))

def death_rate(N, S, d_0, beta = 1, gamma=1):
    return d_0 * (1/2) * ((1 / (1+beta*S)) + (N / (N + gamma*S)))

def death_rate_no_N_depen(S, d_0, beta = 1):
    return d_0 *  ((1 / (1+beta*S)) )



#scales, N_0s = get_random_joint_kde_samples(20)
#divisions = sample_divisions(20)
#print(scales)

#N_0 = round(N_0s[1])
N_0 = 1000000000

#scale = scales[1]
#d_0 = 1/scale
d_0 = 0.004
#division = divisions[1]

c_birth = 1.1#N_0/division#1#N_0 /division * 0.0001
c_death = 0.005
B = 0.01#0.01
S = 0

t_steps = 1000

N_t = []
T_t = []


for x in range(t_steps):
    if N_0 <= 10**3:
        break
    #N_death_0 = np.random.poisson(N_0 * d_0)
    N_deaths = np.random.poisson(N_0 * death_rate_no_N_depen(S=S, d_0=d_0, beta = 10000) )

    #else:
    #    N_deaths = np.random.poisson(N_0 * shape )
    #print(N_0, N_deaths)
    N_0 -= N_deaths
    S += N_deaths * B
    S -= N_deaths * c_death

    print(S)

    #if S > 0:
    #    #N_births = np.random.poisson(N_0 * (birth_rate_no_N_depen(N_0, S, c_birth)))

    #    N_births = np.random.poisson(N_0 * (birth_rate_no_N_depen(S=S, c=c_birth)))
    #    N_0 += N_births
    #    S -= N_births * c_birth

    N_t.append(N_0)
    T_t.append(x)



fig = plt.figure()
plt.scatter(T_t, N_t, color = '#175ac6')
plt.ylim([min(N_t)*0.5,max(N_t)*2])
plt.yscale('log')



plt.xlabel("Time", fontsize = 16)
plt.ylabel("Population size", fontsize = 14)
fig.tight_layout()
fig.savefig(os.path.expanduser("~/GitHub/LTDE") + '/figs/test_death_sim.png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
