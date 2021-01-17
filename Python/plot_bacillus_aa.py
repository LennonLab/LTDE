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



df = pd.read_csv(lt.get_path() + '/data/mzML_Files_forR/Bacillus_AA_Conc_1000d.csv', sep = ',')
df = df.set_index('Name')
aa_list = ['Ala', 'Gly', 'Val', 'Leu', 'Ile', 'Pro', 'Met', 'Ser', 'Thr', 'Phe',
        'Asp', 'Glu', 'Orn', 'Lys', 'His', 'Tyr', 'Cys-Cys']
aa_dict = {'Ala':'Alanine', 'Gly':'Glycine', 'Val':'Valine', 'Leu':'Leucine',
            'Ile': 'Isoleucine', 'Pro':'Proline', 'Met':'Methionine',
            'Ser':'Serine', 'Thr':'Threonine', 'Phe':'Phenylalanine',
            'Asp':'Aspartic Acid', 'Glu':'Glutamic acid', 'Orn':'Arginine',
            'Lys':'Lysine', 'His':'Histidine', 'Tyr':'Tyrosine','Cys-Cys':'Cystine'}
molar_mass_dict = {'Ala':89.094, 'Gly':75.07, 'Val':117.2, 'Leu':113.2,
                    'Ile':113.2, 'Pro':115.1, 'Met':149.2, 'Ser':105.1,
                    'Thr':119.1, 'Phe':165.2, 'Asp':133.1, 'Glu':147.1,
                    'Orn':174.2, 'Lys':146.2, 'His':155.2, 'Tyr':181.2,
                    'Cys-Cys':240.1}

bio_reps = ['KBS0812A', 'KBS0812B', 'KBS0812C', 'KBS0812D']
fig = plt.figure()
plt.axvline(0, color = 'dimgrey', lw = 2, ls = '--', zorder=1)
for i, aa in enumerate(aa_list):
    aa_series = df[df.columns[df.columns.to_series().str.contains(aa)]]
    column_name = aa_series.columns.values[0]
    blank = aa_series.loc[ 'std8' , : ].values[0]
    aa_series[column_name] = aa_series[column_name].apply(lambda x: x - blank)
    for bio_rep in bio_reps:
        aa_series_bio_rep = aa_series.loc[ df.index[df.index.to_series().str.contains(bio_rep)] , : ]
        aa_measures = aa_series_bio_rep[column_name].values
        aa_measures = (aa_measures*0.0001) / molar_mass_dict[aa]
        if len(aa_measures) < 3:
            continue
        mean = np.mean(aa_measures)
        sem = np.std(aa_measures) / math.sqrt(len(aa_measures))

        jitter_rep = np.random.normal(i, 0.08)

        plt.errorbar(mean, i, xerr = sem*2, \
            fmt = 'o', alpha = 0.7, barsabove = True, marker = '.', \
            mfc = 'b', mec = 'none', c = 'k', zorder=3, ms=17)

y1 = list(range(len(aa_list)))
latex_labels = [aa_dict[x] for x in aa_list]

plt.xlabel("Blank-corrected molar concentration, mol/L" , fontsize = 12)

plt.yticks(y1, latex_labels, rotation=0)
fig.savefig(lt.get_path() + '/figs/bacillus_aa.pdf', format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
