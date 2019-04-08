from __future__ import division
import os, math
import pandas as pd
import numpy as np
import  matplotlib.pyplot as plt
import ltde_tools as lt
from sklearn.model_selection import GridSearchCV
from sklearn.neighbors import KernelDensity


def fig1():
    df_path = lt.get_path() + '/data/demography/longtermdormancy_20170620_nocomments.csv'
    df = pd.read_csv(df_path, sep = ',')

    df['N'] = df['Colonies'] * (10 ** (df['Dilution'] + 1))
    df['Dormstart_date'] =  pd.to_datetime(df['Dormstart_date'], format='%d-%b-%y')
    df['Firstread_date'] =  pd.to_datetime(df['Firstread_date'], format='%d-%b-%y')

    df_KBS0714 = df.loc[df['Strain'] == 'KBS0714']
    df_KBS0715 = df.loc[df['Strain'] == 'KBS0715']
    df_KBS0812 = df.loc[df['Strain'] == 'KBS0812']

    day0_KBS0714 = pd.to_datetime({'year':[2013], 'month':[1], 'day':[30]})
    df_KBS0714['Days'] = df_KBS0714['Firstread_date'] - day0_KBS0714[0]
    days_KBS0714 = (df_KBS0714['Days'] / np.timedelta64(1, 'D')).astype(int)

    day0_KBS0715 = pd.to_datetime({'year':[2013], 'month':[2], 'day':[11]})
    df_KBS0715['Days'] = df_KBS0715['Firstread_date'] - day0_KBS0715[0]
    days_KBS0715 = (df_KBS0715['Days'] / np.timedelta64(1, 'D')).astype(int)

    day0_KBS0812 = pd.to_datetime({'year':[2013], 'month':[3], 'day':[8]})
    df_KBS0812['Days'] = df_KBS0812['Firstread_date'] - day0_KBS0812[0]
    days_KBS0812 = (df_KBS0812['Days'] / np.timedelta64(1, 'D')).astype(int)

    df_mean_death = lt.get_mean_time_death()
    df_mean_death_KBS0714 = df_mean_death.loc[df_mean_death['strain'] == 'KBS0714']
    df_mean_death_KBS0715 = df_mean_death.loc[df_mean_death['strain'] == 'KBS0715']
    df_mean_death_KBS0812 = df_mean_death.loc[df_mean_death['strain'] == 'KBS0812']
    df_mean_death_KBS0714_alpha = np.mean(df_mean_death_KBS0714.alpha.values)
    df_mean_death_KBS0715_alpha = np.mean(df_mean_death_KBS0715.alpha.values)
    df_mean_death_KBS0812_alpha = np.mean(df_mean_death_KBS0812.alpha.values)
    df_mean_death_KBS0714_beta = np.mean(df_mean_death_KBS0714.beta.values)
    df_mean_death_KBS0715_beta = np.mean(df_mean_death_KBS0715.beta.values)
    df_mean_death_KBS0812_beta = np.mean(df_mean_death_KBS0812.beta.values)

    days = [days_KBS0714, days_KBS0715, days_KBS0812]
    Ns = [df_KBS0714.N.values, df_KBS0715.N.values, df_KBS0812.N.values]

    fig = plt.figure()

    y_axis_range = [0.5,10000000000]

    weib_x = list(range(0, 1001))

    def weibull_predict(alpha, beta, N0, t):
        return N0 * math.exp(- ((t/beta) ** alpha) )

    weib_y_714 = [weibull_predict(df_mean_death_KBS0714_alpha, df_mean_death_KBS0714_beta, Ns[0][0], x) for x in weib_x]
    weib_y_715 = [weibull_predict(df_mean_death_KBS0715_alpha, df_mean_death_KBS0715_beta, Ns[1][0], x) for x in weib_x]
    weib_y_812 = [weibull_predict(df_mean_death_KBS0812_alpha, df_mean_death_KBS0812_beta, Ns[2][0], x) for x in weib_x]


    ax1 = plt.subplot2grid((3, 3), (0, 0), colspan=1)
    ax1.scatter(days_KBS0714, df_KBS0714.N, alpha=0.1, color = '#1f77b4', zorder=1)
    ax1.plot(weib_x, weib_y_714, c = 'k', zorder=2, linestyle = '--')
    ax1.set_yscale("log")
    ax1.set_xlim([-100,1100])
    ax1.set_ylim(y_axis_range)
    ax1.set_ylabel('CFUs', fontsize = 14)
    ax1.set_title(r"${Micrococcus} \; \mathrm{sp. KBS0714}$", fontsize = 9)
    ax1.tick_params(axis = 'y', which = 'major', labelsize = 6)
    ax1.text(700, 500000000, r"${\alpha}  = " + str(round(df_mean_death_KBS0714_alpha,2)) + "$" , fontsize = 7)
    ax1.text(700, 20000000, r"${\beta}  = " + str(round(df_mean_death_KBS0714_beta,2)) + "$" , fontsize = 7)


    ax2 = plt.subplot2grid((3, 3), (0, 1), colspan=1)
    ax2.scatter(days_KBS0715, df_KBS0715.N, alpha=0.1, color = '#1f77b4', zorder=1)
    ax2.plot(weib_x, weib_y_715, c = 'k', zorder=2, linestyle = '--')
    ax2.set_xlim([-100,1100])
    ax2.set_ylim(y_axis_range)
    ax2.set_yscale("log")
    ax2.set_xlabel('Time (days)', fontsize = 14)
    ax2.set_title(r"${Arthrobacter} \; \mathrm{sp. KBS0715}$", fontsize = 9)
    ax2.tick_params(axis = 'y', which = 'major', labelsize = 6)
    ax2.set_yticklabels([])
    ax2.text(475, 500000000, r"${\alpha}  = " + str(round(df_mean_death_KBS0715_alpha,2)) + "$" , fontsize = 7)
    ax2.text(475, 20000000, r"${\beta}  = " + str(round(np.log10(df_mean_death_KBS0715_beta),2)) + "\, , \mathrm{log_{10}}$" , fontsize = 7)
    #print(df_mean_death_KBS0715_beta)

    ax3 = plt.subplot2grid((3, 3), (0, 2), colspan=1)
    ax3.scatter(days_KBS0812, df_KBS0812.N, alpha=0.1, color = '#1f77b4', zorder=1)
    ax3.plot(weib_x, weib_y_812, c = 'k', zorder=2, linestyle = '--')
    ax3.set_xlim([-100,1100])
    ax3.set_ylim(y_axis_range)
    ax3.set_yscale("log")
    ax3.get_xaxis().get_major_formatter().labelOnlyBase = False
    ax3.minorticks_off()
    ax3.set_title(r"${Bacillus} \; \mathrm{sp. KBS0812}$", fontsize = 9)
    ax3.tick_params(axis = 'y', which = 'major', labelsize = 6)
    ax3.set_yticklabels([])
    ax3.text(550, 500000000, r"${\alpha}  = " + str(round(df_mean_death_KBS0812_alpha,2)) + "$" , fontsize = 7)
    ax3.text(550, 20000000, r"${\beta}  = " + str(round(np.log10(df_mean_death_KBS0812_beta),2)) + "\, , \mathrm{log_{10}}$" , fontsize = 7)


    plt.tight_layout()
    fig_name = lt.get_path() + '/figs/fig1.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()





def fig2():
    df = lt.get_mean_time_death()
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
    fig_name = lt.get_path() + '/figs/fig2.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()




fig1()
#fig2()
