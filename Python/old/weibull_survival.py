from __future__ import division
import os, math
import  matplotlib.pyplot as plt
#import seaborn as sns
import numpy as np
import pandas as pd
from scipy.stats import weibull_min, gamma, poisson
from scipy import stats
from scipy.optimize import minimize
from scipy.misc import factorial


#.weibull_min as wm

mydir = os.path.expanduser("~/GitHub/LTDE/")

def weibull_plot():
    x = np.linspace(0,20,100)
    a = 1.5
    c = 5

    x1 = weibull_min.cdf(x, 0.5, loc=0, scale=c)
    x2 = weibull_min.cdf(x, 1, loc=0, scale=c)
    x3 = weibull_min.cdf(x, 1.5, loc=0, scale=c)


    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(x, x1, 'b:', label=r'$\alpha = 0.5$')
    ax.plot(x, x2, 'k', label=r'$\alpha = 1.0$')
    ax.plot(x, x3, 'r--', label=r'$\alpha = 1.5$')
    plt.legend(loc = 'lower right')

    plt.xlabel(r'$t$', fontsize = 12)
    plt.ylabel(r'$F(t)$', fontsize = 11)

    fig_name = mydir + 'weibull.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()

def weibull_survival_plot():
    x = np.linspace(0,20,100)
    a = 1.5
    c = 5

    x1 = weibull_min.cdf(x, 0.5, loc=0, scale=c)
    x2 = weibull_min.cdf(x, 1, loc=0, scale=c)
    x3 = weibull_min.cdf(x, 1.5, loc=0, scale=c)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(x, np.log(1-x1), 'b:', label=r'$\alpha = 0.5$')
    ax.plot(x, np.log(1-x2), 'k', label=r'$\alpha = 1.0$')
    ax.plot(x, np.log(1-x3), 'r--', label=r'$\alpha = 1.5$')
    plt.legend(loc = 'upper right')

    plt.xlabel('Time (days)', fontsize = 12)
    plt.ylabel('Proportion surviving, ' + r'$ln\,S(t)$', fontsize = 11)

    fig_name = mydir + 'weibull_survival.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()

def half_life_hist():
    df = pd.read_csv(mydir + 'data/weibull_log_results.csv', sep = ',')
    half_life = np.log10(df.Half_life.values)
    fig = plt.figure()
    plt.hist(half_life, normed=True, bins=30)
    plt.xlabel(r'$t_{1/2}$' + ', ' + r'$log_{10}$', fontsize = 12)
    plt.ylabel('Frequency', fontsize = 11)
    fig_name = mydir + 'half_life_hist.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()

def half_life_gamma():
    df = pd.read_csv(mydir + 'data/weibull_log_results.csv', sep = ',')
    half_life = np.log10(df.Half_life.values)
    shape, loc, scale = gamma.fit(half_life)
    x = np.linspace(min(half_life), max(half_life), num =1000)
    print x
    g = gamma.pdf(x=x, a=shape, loc=loc, scale=scale)

    fig = plt.figure()
    plt.hist(half_life, normed=True, bins=30)
    plt.plot(x, g, 'k-', linewidth=6)
    plt.xlabel(r'$t_{1/2}$' + ', ' + r'$log_{10}$', fontsize = 12)
    plt.ylabel('Frequency', fontsize = 11)
    fig_name = mydir + 'half_life_gamma.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()

    #fig, axes = plt.subplots(1, 3, figsize=(13,4))
    #ax = axes[0]
    #ax.hist(y1, bins=40, normed=True);
    #ax.annotate(s='shape = %.3f\nloc = %.3f\nscale = %.3f' %(shape1, loc1, scale1), xy=(6,.2))
    #ax.set_title('gamma fit')




#def half_life_bp():
#    df = pd.read_csv(mydir + 'data/weibull_log_results.csv', sep = ',')
#    df1 = df[['strain','Half_life']]
#    df1['Half_life_log'] = np.log10(df1.Half_life)

#    f, ax = plt.subplots(figsize=(11, 15))

#    #ax.set_axis_bgcolor('#fafafa')
#    #plt.title("Box Plot of Transformed Data Set (Breast Cancer Wisconsin Data Set)")
#    #ax.set(xlim=(-.05, 1.05))
#    #ax = sns.boxplot(data = df1, orient = 'h', palette = 'Set2')
#    ax = sns.boxplot(x="strain", y="Half_life_log", data=df1, orient = "v" )
#    #ax.set(xlabel='', ylabel=r'$t_{1/2},\, log_{10}$')
#    ax.set_xlabel("",fontsize=0)
#    ax.set_ylabel(r'$t_{1/2},\, log_{10}$',fontsize=30)
#    for item in ax.get_xticklabels():
#        item.set_rotation(90)

#    ax.figure.savefig(mydir + 'half_life_bp.png')


def irep_vs_half_life():
    df1 = pd.read_csv(mydir + 'data/weibull_log_results.csv', sep = ',')
    df1['Half_life_log'] = np.log10(df1.Half_life)

    df2 = pd.read_csv(mydir + 'iRepFinal.txt', sep = ' ')
    df2 = df2.drop(df2.index[7])
    df2 = df2.dropna(subset=['iRep'])
    df2.loc[40,'Strain'] = 'KBS0711W'
    df2.loc[41,'Strain'] = 'KBS0711W'
    df2.loc[42,'Strain'] = 'KBS0711W'
    df2.loc[43,'Strain'] = 'KBS0711W'
    df2.loc[40,'Replicate'] = 'A'
    df2.loc[41,'Replicate'] = 'B'
    df2.loc[42,'Replicate'] = 'C'
    df2.loc[43,'Replicate'] = 'D'

    alpha_to_num = {'A':1, 'B':2, 'C':3, 'C1':3, 'D':4, 'E':5, 'F':6, 'K':10, 'L':11}

    def num_rep(row):
        return alpha_to_num[row['Replicate']]

    df2['rep'] = df2.apply(num_rep, axis=1)
    df2.rename(columns={'Strain':'strain'}, inplace=True)

    df_merged = pd.merge(df1, df2, on=['strain', 'rep'])
    x = np.log10(df_merged.a.values)
    y = np.log10(df_merged.iRep.values)

    fig = plt.figure()
    plt.scatter(x, y, c='#87CEEB', marker='o')
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    predict_y = intercept + slope * x
    pred_error = y - predict_y
    degrees_of_freedom = len(x) - 2
    residual_std_error = np.sqrt(np.sum(pred_error**2) / degrees_of_freedom)
    plt.plot(x, predict_y, 'k-', lw = 5, c = 'black', label = '_nolegend_' )
    #plt.plot(x_M, predict_y_M, 'k-', lw = 2, c = '#87CEEB', label = '_nolegend_')

    #plt.xlabel('Mid-parent height (inches)', fontsize = 18)
    #plt.ylabel('Offspring height (inches)', fontsize = 18)
    #plt.legend(loc = 'upper left')
    fig.tight_layout()
    fig_name = mydir + 'irep_vs_half_life.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()


def n_vs_t_fix_mut():
    N = 1000000
    U = 0.00000001
    s = 0.0001
    x = np.logspace(1, 6, num = 20)

    #print N * np.log(N) * U

    t_fix = (2*np.log(x)) / s
    t_mut = 1 / (2*x*U*s)

    fig = plt.figure()

    plt.plot(x, t_fix, 'b:', label=r'$t_{fix}$')
    plt.plot(x, t_mut, 'r--', label=r'$t_{mut}$')

    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel(r'$N$' + ', ' + r'$log_{10}$', fontsize = 14)
    plt.ylabel('Generations, ' + r'$log_{10}$', fontsize = 14)
    plt.legend(loc = 'upper right', fontsize = 14)
    fig.tight_layout()
    fig_name = mydir + 'n_vs_t_fix_mut.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()


def u_vs_t_fix_mut():
    N = 1000000
    U = 0.00000001
    s = 0.0001
    x = np.logspace(-8, -3, num = 20)

    t_fix_i = (2*np.log(N)) / s
    t_fix = np.asarray([t_fix_i] * 20)
    t_mut = 1 / (2*N*x*s)

    fig = plt.figure()

    plt.plot(x, t_fix, 'b:', label=r'$t_{fix}$')
    plt.plot(x, t_mut, 'r--', label=r'$t_{mut}$')

    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel(r'$U$' + ', ' + r'$log_{10}$', fontsize = 14)
    plt.ylabel('Generations, ' + r'$log_{10}$', fontsize = 14)
    plt.legend(loc = 'upper right', fontsize = 14)
    fig.tight_layout()
    fig_name = mydir + 'u_vs_t_fix_mut.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()


def half_life_vs_null():
    df = pd.read_csv(mydir + 'data/weibull_log_results.csv', sep = ',')
    y = np.log10(df.Half_life.values)
    #null_log = np.log10(df.a.values) + ( (1/df.b.values) * np.log10(np.log(0.5) * -1) )
    x = np.log10(df.a.values) +  np.log10(np.log(0.5) * -1)
    fig = plt.figure()

    plt.scatter(x, y, c='#87CEEB', marker='o')
    plt.plot([-5, 11],[-5, 11], '--', color = 'dimgrey')

    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    predict_y = intercept + slope * x
    pred_error = y - predict_y
    degrees_of_freedom = len(x) - 2
    residual_std_error = np.sqrt(np.sum(pred_error**2) / degrees_of_freedom)
    plt.plot(x, predict_y, 'k-', lw = 1.5, c = 'black', label = 'm = ' + str(slope) )

    print slope
    print r_value ** 2

    plt.xlabel('Half-life with linear decay\n' + r'$\left ( \alpha=  1 \right )$', fontsize = 14)
    plt.ylabel('Half-life allowing for non-linear decay\n' + r'$\left ( \alpha\neq  1 \right )$', fontsize = 14)

    fig.tight_layout()
    fig_name = mydir + 'output/half_life_vs_null.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()

    #half_life_null = np.log10(df.Beta.values)





#half_life_vs_null()
half_life_hist()
