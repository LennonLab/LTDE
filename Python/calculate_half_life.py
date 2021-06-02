from __future__ import division

import glob, math, os, re
import pandas as pd
import numpy as np
import ltde_tools as lt
from scipy import stats
from scipy.stats import t
from scipy.integrate import odeint
from decimal import Decimal
import _pickle as pickle

#from sklearn.model_selection import GridSearchCV
#from sklearn.neighbors import KernelDensity
import statsmodels.stats.multitest as mt
import statsmodels.formula.api as smf

from Bio import SeqIO

from statsmodels.base.model import GenericLikelihoodModel


filepath = lt.get_path() + '/data/demography/weibull_results_clean.csv'

half_life_dict = {}

count = 0
for line in open(filepath, 'r'):

    if count == 0:
        count += 1
        continue

    line_split = line.strip().split(',')

    taxon = line_split[1]

    #print(line_split[3])

    d_0 = float( line_split[3])
    k = float( line_split[4])



    half_life = (1/d_0) * ((-1* np.log(0.5))**(1/k) )
    #T_{ext} =\ d_{0}^{-1}(-\mathrm{log_{e}}{(S_{ext})}^{1/k})


    half_life_years = half_life /365

    if taxon not in half_life_dict:
        half_life_dict[taxon] = []


    half_life_dict[taxon].append(half_life_years)

    count += 1


for key in half_life_dict.keys():

    mean_half_life = 10**np.mean(np.log10(half_life_dict[key]))


    print(key, np.log10(mean_half_life))



# https://doi.org/10.1093/aob/mcp082
# days
plant_half_lives = [76.9,38.7,54.4,52.6,47.0,43.6,11.2,8.3,39.8,0.1,78.2,69.3,57.0,54.0,10.7,47.0,59.3,9.4,83.5,19.2,82.8,0.9,41.7,49.5,53.7,59.2,5.6,69.0,9.0,14.0,4.5,53.9,36.7,26.7,22.0,38.7,51.1,31.0,14.9,21.8,16.0,4.6,25.4,62.9,84.4,6.3]

plant_half_lives = np.asarray(plant_half_lives)

plant_half_lives  = plant_half_lives/ 365

#print()


print(max(np.log10(plant_half_lives)))
