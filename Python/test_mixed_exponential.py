from __future__ import division

import glob, math, os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import ltde_tools as lt


N1 = 10**7
N2 = 10**7


def calculate_bi_exponential(N1, N2, d1, d2, t):

    return (N1*np.exp(-1*d1*t)) + (N2*np.exp(-1*d2*t))



d1_list = [0.00001, 0.0001, 0.001, 0.01, 0.1]
d2 = 0.001

fig = plt.figure()

colors = ['red', 'tomato', 'grey', 'lightblue', 'dodgerblue']

for d1_idx, d1 in enumerate(d1_list):

    N_list = []

    t_list = list(range(1000))

    for t in t_list:

        N_list.append(calculate_bi_exponential(N1, N2, d1, d2, t))


    print(colors[d1_idx])

    plt.plot(t_list, N_list, zorder=2,ls='--', label= r'$d_{1}/d_{2}=$' + str(d1/d2), c=colors[d1_idx], lw=2)



plt.xlabel('Days, ' + r'$t$', fontsize = 16)
plt.ylabel('Population size, ' + '$N(t)$', fontsize = 16)
plt.yscale('log', base=10)


plt.legend(loc='upper right', prop={'size': 8})

#fig.savefig(lt.get_path() + '/figs/spoiie_death_curve.pdf', format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
fig.savefig(lt.get_path() + '/figs/test_exponential.pdf', format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)

plt.close()
