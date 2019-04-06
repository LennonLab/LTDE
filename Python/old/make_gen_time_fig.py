from __future__ import division
import os, math
import pandas as pd
import numpy as np
import  matplotlib.pyplot as plt
import ltde_tools as lt





def get_weighted_mean_time_death():
    df = lt.get_mean_time_death()
    df['strain'] = df['strain'].str.replace('KBS0711W','KBS0711')
    # get rid of the insanely out there bacillus estimate
    df = df.drop(df.index[[96, 99]])
    df_strain = df.groupby(['strain'])
    df_out = open(lt.get_path() + '/data/demography/mean_weibull_results.txt', 'w')
    header = ['Strain', 'Number_reps', 'Weighted_mean', 'CI025_weighted_mean', 'CI975_weighted_mean']
    df_out.write('\t'.join(header) + '\n')
    for name, group in df_strain:
        n = group['N.obs'].values
        weighted_mean = sum(group['N.obs'].values * group.mean_days_death.values) /  sum(group['N.obs'].values)
        pooled_sd =  math.sqrt(sum((group['N.obs'].values -1) * (group.sd_days_death.values ** 2)) / sum(group['N.obs'].values -1))
        #print(name, weighted_mean, pooled_sd)
        CI025_weighted_mean_days_death = lt.weibull_CIs(mean = weighted_mean, sd = pooled_sd, n =n, lower = True, pooled = True)
        CI975_weighted_mean_days_death = lt.weibull_CIs(mean = weighted_mean, sd = pooled_sd, n =n, lower = False, pooled = True)
        line = [name, str(group.shape[0]), str(weighted_mean), str(CI025_weighted_mean_days_death), str(CI975_weighted_mean_days_death)]
        df_out.write('\t'.join(line) + '\n')

    df_out.close()


def plot_weighted_mean():
    # merge trait dataset with mean death time
    df_path = lt.get_path() + '/data/demography/mean_weibull_results.txt'
    df = pd.read_csv(df_path, sep = '\t', index_col = 0)
    df_traits_path = lt.get_path() + '/data/traits/persistence.phylo.txt'
    df_traits = pd.read_csv(df_traits_path, sep = '\t', index_col = 0)
    df_traits = df_traits[np.isfinite(df_traits['umax'])]
    df_traits.rename(columns={'Code':'Strain'}, inplace=True)
    # double check with jay that umax is in days
    df_traits['division_time_days'] = df_traits.apply(lambda row: math.log(2)/ row.umax, axis=1)
    df_merge = pd.merge(df_traits, df, left_index=True, right_index=True)
    df_merge['division_time_starv_mean'] = df_merge['division_time_days'] + df_merge['Weighted_mean']
    df_merge['division_time_starv_CI025'] = df_merge['division_time_days'] + df_merge['CI025_weighted_mean']
    df_merge['division_time_starv_CI975'] = df_merge['division_time_days'] + df_merge['CI975_weighted_mean']
    df_merge.sort_values("division_time_starv_mean", inplace=True)
    fig = plt.figure()
    count = 0
    strains = []
    for key, row in df_merge.iterrows():
        if (key == 'KBS0715') or (key == 'KBS0713'):
            continue
        #div_time = np.log10(row['division_time_days'] * 365 * (3.8 * (10**12) ))
        #div_time_starv_CI025 = np.log10(row['division_time_starv_CI025']* 365 * (3.8 * (10**12) ))
        #div_time_starv_CI975 = np.log10(row['division_time_starv_CI975']* 365 * (3.8 * (10**12) ))
        #div_time_starv_mean = np.log10(row['division_time_starv_mean']* 365 * (3.8 * (10**12) ))

        div_time = np.log10(   ( 1 / (row['division_time_days'] / 365)) * (3.8 * (10**12) ))
        div_time_starv_CI025 = np.log10( (1 / (row['division_time_starv_CI025']/ 365)) * (3.8 * (10**12) ))
        div_time_starv_CI975 = np.log10( (1 / (row['division_time_starv_CI975']/ 365)) * (3.8 * (10**12) ))
        div_time_starv_mean = np.log10(  (1 / (row['division_time_starv_mean'] / 365)) * (3.8 * (10**12) ))

        plt.hlines(count, div_time, div_time_starv_CI025, linestyles = ':')
        plt.hlines(count, div_time_starv_CI025, div_time_starv_CI975, linestyles = '-')
        plt.scatter(div_time_starv_mean, count, s=30, c = 'k')
        plt.scatter(div_time, count, s=30, facecolors='none', edgecolors='k')
        count += 1
        strains.append(key)

    plt.yticks(np.arange(0, count), strains)
    plt.vlines(12, 0, count, linestyles = '--')
    plt.xlabel('Number of generations, ' + r'$log_{10}$', fontsize = 14)
    fig.tight_layout()
    fig_name = lt.get_path() + '/figs/range_division_time.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()

    #df_emp_path = lt.get_path() + '/data/EMP/EMP-clean/estimated_cell_counts.txt'
    #df_emp = pd.read_csv(df_emp_path, sep = '\t', index_col = 0)
    #print(df_emp)

    # figure out how to calculate confidence intervals




get_weighted_mean_time_death()
#plot_weighted_mean()
