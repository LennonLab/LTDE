
import os
import pwlf
import numpy as np
import pandas as pd
import ltde_tools as lt

df_colors = pd.read_csv(lt.get_path() + '/data/colors.csv', sep = ',')


def piecewise_regression():
    df = pd.read_csv(os.path.expanduser("~/GitHub/LTDE") + '/data/demography/longtermdormancy_20190528_nocomments.csv', sep = ',')
    # KBS0721 rep
    df['N'] = (df['Colonies']+1) * (1000 / df['Inoculum']) * (10 ** (df['Dilution'] ))

    df['Dormstart_date'] =  pd.to_datetime(df['Dormstart_date'], format='%d-%b-%y')
    df['Firstread_date'] =  pd.to_datetime(df['Firstread_date'], format='%d-%b-%y')
    df['Days'] = df['Firstread_date'].sub(df['Dormstart_date'], axis=0)
    df['Days'] = df['Days'].dt.days.astype('int')
    #sdf
    N_dead_list = []
    delta_slope_list = []
    time_split = []
    slope1_list = []
    slope2_list = []
    taxa = list(set(df.Strain.to_list()))
    #fig_all = plt.figure()
    flux_list = []
    slope2_scale_list = []

    df_out = open(lt.get_path() + '/data/demography/piecewise_regression.txt', 'w')
    df_out.write('\t'.join(['Species', 'rep', 'N0', 'slope1' , 'slope2', 'time_split', 'N_split']) + '\n')

    for taxon in taxa:
        print(taxon)
        if taxon == 'KBS0714':
            continue
        if taxon == 'KBS0719':
            continue
        if taxon == 'KBS0711W':
            continue
        if taxon == 'KBS0816':
            continue
        if taxon == 'KBS0727':
            continue
        if taxon == 'KBS0704':
            continue

        #if taxon == 'KBS0725':
        #    continue
        #if taxon == 'KBS0702':
        #    continue
        #if taxon == 'KBS0712':
        #    continue
        #if taxon == 'KBS0703':
        #    continue
        #if taxon == 'KBS0706':
        #    continue

        taxon_color = df_colors.loc[df_colors['strain'] == taxon].Color.to_list()[0]

        #if taxon != 'KBS0710':
        #    continue
        #fig = plt.figure()
        df_taxon = df[(df["Strain"] == taxon)]
        reps = list(set(df_taxon.Rep.to_list()))
        for rep in reps:
            df_taxon_rep =  df_taxon[(df_taxon["Rep"] == rep)]
            df_taxon_rep.sort_values('Days')

            x = df_taxon_rep.Days.values
            if len(x) < 20:
                continue
            y = np.log10(df_taxon_rep.N.values)

            N0 = df_taxon_rep.N.values[0]

            my_pwlf = pwlf.PiecewiseLinFit(x, y)


            # fit the data for four line segments
            res = my_pwlf.fit(2)

            # predict for the determined points
            xHat = np.linspace(min(x), max(x), num=10000)
            yHat = my_pwlf.predict(xHat)

            #print(taxon, len(res))

            N_switch = my_pwlf.intercepts[0] + (my_pwlf.calc_slopes()[0] *res[1]  )

            #N_flux =

            slopes = my_pwlf.calc_slopes()
            N_dead = (10 ** max(y)) - (10 ** N_switch)
            N_dead_flux = ((10 ** max(y)) * slopes[0])
            #angle_switch = np.arctan(np.absolute( (slopes[1]-slopes[0]) /(1+ (slopes[0]*slopes[1])) ))
            delta_slope = (slopes[1] - slopes[0])
            #print(taxon, rep, N_switch, delta_slope)

            N_dead_list.append(N_dead / res[1])
            #N_dead_per_time_list.append(N_dead / )
            slope1_list.append(slopes[0])
            slope2_list.append(slopes[1])
            delta_slope_list.append(delta_slope)

            time_split.append(res[1])

            #taxon_color
            #print(delta_slope)
            #max(10**y) * slopes[0]





            #plt.scatter(np.log10(np.abs(N_dead_flux / res[1] ) ), np.log10(delta_slope), c = taxon_color)
            #plt.scatter(np.log10(np.abs(slopes[0]) ), np.log10(slopes[1]), c = taxon_color)
            flux = np.abs( (10**my_pwlf.intercepts[0]) * slopes[0])
            flux_list.append(flux)

            #if taxon == 'KBS0812':
            #    print(flux, my_pwlf.intercepts[0], slopes[1])

            slope2_scale_list.append(slopes[1] - slopes[0])


            #df_out.write('\t'.join([ strain, str(round(weighted_mean_cov, 3)) ]) + '\n')


            df_out.write('\t'.join([taxon, str(rep), str(N0), str(slopes[0]) , str(slopes[1]), str(res[1]), str(10**N_switch)]) + '\n')


            #plt.scatter(np.log10(flux ), np.log10(slopes[1] ), c = taxon_color)

            #plt.scatter(slopes[0], slopes[1], c = taxon_color)


            #plt.scatter(x, y)
            #plt.plot(xHat, yHat, '-')


        #plt.xscale('log',basex=10)
        #plt.yscale('log',basey=10)
        #plt.xlim(10**-2,10**10)
        #plt.xlabel('Number of dead cells', fontsize = 12)
        #plt.ylabel('absolute value of second slope, log10', fontsize = 12)
        #fig.savefig(lt.get_path() + '/figs/taxon_piece/'+taxon+'.png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
        #plt.close()

    df_out.close()


    #plt.xscale('log', basex=10)
    #plt.yscale('log', basey=10)
    #plt.yscale('symlog')

    #plt.xlim(-0.1,0)

    #print(np.log10(flux_list), np.log10(slope2_scale_list))
    #slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(flux_list), slope2_scale_list)

    #print(slope, p_value)

    #plt.scatter(np.log10(N_dead_list), np.log10(np.abs(delta_slope_list)))
    #plt.xlabel('Flux N(0)*Slope1 , log10', fontsize = 12)
    #plt.ylabel('Slope2', fontsize = 12)
    #plt.plot(xHat, yHat, '-')
    #fig_all.savefig(lt.get_path() + '/figs/piece.png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    #plt.close()


piecewise_regression()
