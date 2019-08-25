from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import os

N = 10000000
d_0 = 0.001
S_0 = 0
alpha = 100
N_t = []
S_min = 0.001
alpha = 0.0001
# maintenance min = 0.01 × 10−19 J s−1 cell−1
# power from necromass oxidation
t_steps = 10000
for x in range(t_steps):
    print(N)

    d_t = d_0 * np.exp(-1 * max(0, (S_0/N)-S_min))
    #print( d_t)
    delta_N = np.random.poisson(d_t*N)
    N -= delta_N
    S_0 += delta_N * alpha
    N_t.append(N)


fig = plt.figure()
plt.scatter(list(range(t_steps)), N_t, color = '#175ac6')
plt.yscale('log')

plt.xlabel("Time", fontsize = 16)
plt.ylabel("Population size", fontsize = 14)
fig.tight_layout()
fig.savefig(os.path.expanduser("~/GitHub/LTDE") + '/figs/test_death_sim.png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
