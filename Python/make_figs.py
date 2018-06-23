from __future__ import division
import os, math
import pandas as pd
import numpy as np
import  matplotlib.pyplot as plt
from matplotlib import rc

import ltde_tools as lt
from sklearn.model_selection import GridSearchCV
from sklearn.neighbors import KernelDensity
from operator import itemgetter


def fig1():
    df = lt.get_mean_time_death()
    df = lt.clean_demography_df(df)
    fig = plt.figure()
    # alpha kde
    alpha = df.alpha.values
    grid_alpha = GridSearchCV(KernelDensity(),
                    {'bandwidth': np.linspace(0.1, 10, 50)},
                    cv=20) # 20-fold cross-validation
    grid_alpha.fit(alpha[:, None])
    x_grid_alpha = np.linspace(0, 2.5, 1000)
    kde_alpha = grid_alpha.best_estimator_
    pdf_alpha = np.exp(kde_alpha.score_samples(x_grid_alpha[:, None]))
    ax1 = plt.subplot2grid((2, 2), (0, 0), colspan=1)
    pdf_alpha = [x / sum(pdf_alpha) for x in pdf_alpha]
    ax1.plot(x_grid_alpha, pdf_alpha, alpha=0.8, lw = 2, color = '#1f77b4') #, marker='o')
    #ax1.fill_between(x_grid_alpha, 0, pdf_alpha)
    #ax1.hist(alpha, 35, fc='gray', histtype='stepfilled', alpha=0.3, normed=True)
    ax1.axvline(x=1, color='darkgrey', linestyle='--', lw = 2.5)
    ax1.axvline(x=np.mean(alpha), color='#1f77b4', linestyle='--', lw = 2.5)
    ax1.set_xlim([-0.1,2.6])
    ax1.set_xlabel("Scale parameter, " +r'$\alpha$', fontsize = 14)
    ax1.set_ylabel("Probability density", fontsize = 14)

    # half life kde
    half_life = np.log10(df.half_life.values)
    grid_half_life = GridSearchCV(KernelDensity(),
                    {'bandwidth': np.linspace(0.1, 10, 50)},
                    cv=20) # 20-fold cross-validation
    grid_half_life.fit(half_life[:, None])
    x_grid_half_life = np.linspace(-10, 15, 1000)
    kde_half_life = grid_half_life.best_estimator_
    pdf_half_life = np.exp(kde_half_life.score_samples(x_grid_half_life[:, None]))
    pdf_half_life = [x / sum(pdf_half_life) for x in pdf_half_life]
    ax2 = plt.subplot2grid((2, 2), (1, 0), colspan=1)
    ax2.plot( x_grid_half_life, pdf_half_life , color="orange", alpha=0.8)
    ax2.axvline(x=np.mean(half_life), color='orange', linestyle='--', lw = 2.5)
    ax2.set_xlim([-11,11])
    ax2.set_xlabel("Half-life " + r'$\mathrm{(d^{-1}), \, log_{10}} $', fontsize = 14)
    ax2.set_ylabel("Probability density", fontsize = 14)

    ax3 = plt.subplot2grid((2, 2), (0, 1), rowspan=2)
    strains = list(set(df.strain.values))
    half_lives = [np.log10(df.loc[df['strain'] == strain].half_life.values) for strain in strains]
    mean_half_life = [np.median(np.log10(df.loc[df['strain'] == strain].half_life.values)) for strain in strains]

    zipped_half_lives = list(zip(strains, half_lives, mean_half_life))
    zipped_half_lives.sort(key=lambda x: int(x[2])) # reverse= True)
    zipped_half_lives_sorted = sorted(zipped_half_lives, key=lambda x: x[2])
    strains = [x[0] for x in zipped_half_lives_sorted]
    half_lives = [x[1] for x in zipped_half_lives_sorted]
    mean_half_life = [x[2] for x in zipped_half_lives_sorted]

    strain_dict = lt.get_strain_genus_dict()
    genera_labels = list(reversed([strain_dict[x] for x in strains]))
    strain_labels = list(reversed([' sp. ' + x for x in strains]))
    #print(strain_labels)

    #genera
    ax3.boxplot(half_lives, vert = False)
    ax3.yaxis.set_major_formatter(plt.NullFormatter())
    ax3.set_xlim([-6,11])
    #ax3.yaxis.set_ticks_position('none')
    #ax3.gca().xaxis.set_major_locator(plt.NullLocator())
    #ax3.axis('off')
    ax3.set_xlabel("Half-life " + r'$\mathrm{(d^{-1}), \, log_{10}} $', fontsize = 14)
    for i in range(len(genera_labels)):
        genera_label = genera_labels[i]
        strain_label = strain_labels[i]
        y = len(genera_labels) - i - 0.1
        if i == 0:
            ax3.text(-5, y, r"${" + genera_label + "} \, \mathrm{" + strain_label + "}$", fontsize = 5.5)
        else:
            if genera_label == 'Janthinobacterium':
                fontsize = 5.2
            else:
                fontsize = 5.5
            ax3.text(3.2,  y, r"${" + genera_label + "} \, \mathrm{" + strain_label + "}$", fontsize = fontsize)

    plt.tight_layout()
    fig_name = lt.get_path() + '/figs/fig1.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()




fig1()
