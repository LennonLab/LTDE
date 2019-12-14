from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import os

N = 1.38 * 10**9
d_0 = 0.001
S = 0
alpha = 10
#S_min = 0.0001
S_min = 14.4 * 10**-10
#alpha = 0.0001
alpha = 245 * 10 ** -13
#beta = 0.0001
# maintenance min = 0.01 × 10^−19 J s−1 cell−1
# =  0.01 × 10^−10 nJ s−1 cell−1
# * 24 * 60
# = 14.4 × 10^−10 nJ d−1 cell−1

# power from necromass oxidation
#35 nJ cm^-3 yr^-1,
# times cell volume
# 0.7 μm^3 = 7 *10^-13 cm^3
# = 245 x 10^-13 nJ cell^-1 yr^-1


t_steps = 1000

N_range = [int(float('1e7')), int(float('1e9'))]
lamda_range = [1/100, 1/0.001]

N_t = []
T_t = []

for x in range(t_steps):
    if N <= 0:
        break
    print(N)
    #print(N)
    #d_t = d_0 * np.exp(-1 * max(0, (S-S_min)/N))
    d_t = d_0 * np.exp(-1 *  S/(N*S_min) )
    #print(S)
    delta_N = np.random.poisson(d_t*N)
    N -= delta_N
    #S += max(((delta_N * alpha) - (N*S_min)), 0)
    S += (delta_N * alpha) -  min((N*S_min),S )

    N_t.append(N)
    T_t.append(x)


fig = plt.figure()
plt.scatter(T_t, N_t, color = '#175ac6')
plt.ylim([1, max(N_t)])

plt.yscale('log')




plt.xlabel("Time", fontsize = 16)
plt.ylabel("Population size", fontsize = 14)
fig.tight_layout()
fig.savefig(os.path.expanduser("~/GitHub/LTDE") + '/figs/test_death_sim.png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
