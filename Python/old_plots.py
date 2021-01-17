



def plot_weib_indiv_taxon():
    to_remove_KBS0711 = [10,11,12]
    all_taxa = ['KBS0703', 'ATCC13985', 'ATCC43928', 'KBS0701', 'KBS0702',
                'KBS0705', 'KBS0706', 'KBS0707', 'KBS0710', 'KBS0711',
                'KBS0712', 'KBS0713', 'KBS0714', 'KBS0715', 'KBS0721',
                'KBS0722', 'KBS0724', 'KBS0725', 'KBS0801', 'KBS0802',
                'KBS0812']
    df_counts = pd.read_csv(lt.get_path() + '/data/demography/longtermdormancy_20190528_nocomments.csv', sep = ',')
    df_counts['Abund'] = (df_counts.Colonies.values+1) * (1000 / df_counts.Inoculum.values) * ( 10 ** df_counts.Dilution.values )
    df_counts['Dormstart_date'] = pd.to_datetime(df_counts['Dormstart_date'])
    df_counts['Firstread_date'] = pd.to_datetime(df_counts['Firstread_date'])
    df_counts['Days'] = df_counts['Firstread_date'] - df_counts['Dormstart_date'] + dt.timedelta(days=1)
    df_stats = pd.read_csv(lt.get_path() + '/data/demography/weibull_results_clean.csv', sep = ',')
    for taxon in all_taxa:
        print(taxon)
        #taxon = 'KBS0812'
        fig = plt.figure()
        fig.subplots_adjust(hspace=0.35, wspace=0.35)
        taxon_color = df_colors.loc[df_colors['strain'] == taxon].Color.to_list()[0]
        taxon_genus = df_colors.loc[df_colors['strain'] == taxon].Genus.to_list()[0]
        df_counts_taxon = df_counts.loc[(df_counts['Strain'] == taxon)]
        # ignore 10,11,12
        reps = list(set(df_counts_taxon.Rep.to_list()))
        if taxon == 'KBS0711':
            reps = [j for j in reps if j not in to_remove_KBS0711]
        reps.sort()
        plot_dim = len(reps)/2
        if len(reps) == 4:
            plot_dim_x = 2
            plot_dim_y = 2
            axis_text_font = 8
        elif len(reps) == 6:
            plot_dim_x = 2
            plot_dim_y = 3
            axis_text_font = 6
        else:
            print('Number of reps not recognized')

        for rep in reps:
            df_counts_taxon_rep = df_counts_taxon.loc[(df_counts_taxon['Rep'] == rep)]
            df_stats_taxon_rep = df_stats.loc[(df_stats['rep'] == rep) & (df_stats['strain'] == taxon)]
            scale = df_stats_taxon_rep.beta.to_list()[0]
            shape = df_stats_taxon_rep.alpha.to_list()[0]
            N_0 = df_stats_taxon_rep.N_0.to_list()[0]
            last_day = max(df_counts_taxon_rep.Days.dt.days.to_list())
            time = list(range(0, max(df_counts_taxon_rep.Days.dt.days), 1))
            exp_pred = [ N_0* math.exp(-1* (t / scale) ) for t in time]
            weib_pred = [ N_0* (math.exp(-1* (t / scale)** shape ) )  for t in time]
            ax = fig.add_subplot(plot_dim_x, plot_dim_y, rep)
            ax.scatter(df_counts_taxon_rep.Days.dt.days, df_counts_taxon_rep.Abund.values, c=taxon_color, marker = 'o', s = 70, \
                linewidth = 0.6, alpha = 0.5, zorder=1, edgecolors='none')
            ax.plot(time, exp_pred, zorder=2, c='darkgrey',ls='--', lw=2)
            ax.plot(time, weib_pred, zorder=3, c='k',ls='-', lw=2)
            ax.set_ylim([0.2*min(df_counts_taxon_rep.Abund.values), 4*N_0])

            if (scale > 100) or (scale <0.01):
                scale_plot = "{:.2e}".format(scale)
            else:
                scale_plot = round(scale, 2)

            #ax.text(last_day*0.6, 10**(np.log10(4*N_0 - 0.2*min(df_counts_taxon_rep.Abund.values))*0.9), 'Replicate ' + str(rep), fontsize=axis_text_font, transform=ax.transAxes)
            ax.text(0.6, 0.9, 'Replicate ' + str(rep), fontsize=axis_text_font, transform=ax.transAxes)
            ax.text(0.6, 0.77,  r'$d_{0}=$' + str(1/scale_plot), fontsize=axis_text_font, transform=ax.transAxes)
            ax.text(0.6, 0.64,  r'$k=$' + str(round(shape, 2)), fontsize=axis_text_font, transform=ax.transAxes)

            #ax.text(last_day*0.6, 10**(np.log10(4*N_0 - 0.2*min(df_counts_taxon_rep.Abund.values))*0.83),  r'$\lambda=$' + str(scale_plot), fontsize=axis_text_font, transform=ax.transAxes)
            #ax.text(last_day*0.6, 10**(np.log10(4*N_0 - 0.2*min(df_counts_taxon_rep.Abund.values))*0.76),  r'$k=$' + str(round(shape, 2)), fontsize=axis_text_font, transform=ax.transAxes)

            ax.set_yscale('log')

        fig.suptitle(latex_dict[taxon], fontsize=16)

        fig.text(0.5, 0.02, 'Days, ' + r'$t$', ha='center', fontsize=16)
        fig.text(0.02, 0.5, 'Population size, ' + '$N(t)$', va='center', rotation='vertical', fontsize=16)

        #fig.savefig(lt.get_path() + '/figs/taxon_weibull_100/'+taxon+'.png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
        fig.savefig(lt.get_path() + '/figs/taxon_weibull/'+taxon+'.pdf', format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
        #plt.savefig('destination_path.eps', format='eps')

        plt.close()




def plot_logpvalue_survival():
    fig = plt.figure()
    fig.subplots_adjust(hspace=0.35, wspace=0.35)
    pstar_dict = pickle.load(open(lt.get_path() + '/data/breseq/p_star.txt', 'rb'))
    for i in range(0, len(taxa_to_plot)):
        taxon = taxa_to_plot[i]
        pstar_i = pstar_dict[taxon][1]
        num_significant_i = pstar_dict[taxon][0] -1
        df_path = lt.get_path() + '/data/breseq/logpvalues/' + taxon + '.txt'
        df = pd.read_csv(df_path, sep = '\t', index_col=0)
        new_x = df.P_value.tolist()
        new_obs_y = df.Obs_num.tolist()
        new_null_y = df.Null_num.tolist()

        ax = fig.add_subplot(3, 3, i+1)

        ax.plot(new_x, new_null_y, '-', c='dimgrey', lw=4, alpha = 0.8, zorder=0)
        ax.plot(new_x, new_obs_y, '-', c='royalblue', lw=4, alpha = 0.8, zorder=1)
        if pstar_i <0:
            y_range = [f[1] for f in list(zip(new_x, new_obs_y)) if f[0] > 0]
            ax.plot([1, 1],[5e-02,max(y_range)],'k-',linewidth=0.5, zorder=2)
            ax.plot([-3,1],[max(y_range), max(y_range)],'k-',linewidth=0.5, zorder=3)
            ax.plot([1], [max(y_range)], c='r', marker='o', zorder=4)
        else:
            ax.plot([pstar_i, pstar_i],[5e-02,num_significant_i],'k-',linewidth=0.5, zorder=2)
            ax.plot([-3,pstar_i],[num_significant_i, num_significant_i],'k-',linewidth=0.5, zorder=3)
            ax.plot([pstar_i], [num_significant_i], c='r', marker='o', zorder=4)

        ax.set_xlim([0.25, max(new_x)+1])

        ax.title.set_text(latex_dict[taxon])
        ax.title.set_fontsize(7.5)

        ax.xaxis.set_tick_params(labelsize=6)
        ax.yaxis.set_tick_params(labelsize=6)
        #pvalue_axis.step(observed_ps/log(10), null_pvalue_survival(observed_ps),'-',label='Expected',color='k')
        #pvalue_axis.step(observed_ps/log(10), observed_pvalue_survival,'b-',label='Observed')

    fig.text(0.5, 0.02, '$-\mathrm{log}_{10}P$', ha='center', fontsize=16)
    fig.text(0.02, 0.5, 'Number of genes', va='center', rotation='vertical', fontsize=16)

    fig_name = lt.get_path() + '/figs/logpvalue_survival.pdf'
    fig.savefig(fig_name, format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()







def distance_decay():
    fig = plt.figure()
    fig.subplots_adjust(hspace=0.35, wspace=0.35)
    for i in range(0, len(taxa_to_plot)):
        taxon = taxa_to_plot[i]
        df_mult_path = lt.get_path() + '/data/breseq/mult_genes/' + taxon + '.txt'
        df_mult = pd.read_csv(df_mult_path, sep = ', ', index_col=0)
        gene_names = df_mult.index.to_list()
        mult_dist = {}
        # get gbff
        # GCF_005937985_2_ASM593798v2_genomic.gbff
        dist_dict = {}
        if taxon == 'KBS0712':
            df_gbff_path = lt.get_path() + '/data/genomes/genomes_ncbi/KBS0712/GCF_006323185_2_ASM632318v2_genomic.gbff'
        elif taxon == 'KBS0711':
            df_gbff_path = lt.get_path() + '/data/genomes/genomes_ncbi/KBS0711/GCF_005937955_2_ASM593795v2_genomic.gbff'
        elif taxon == 'ATCC13985':
            df_gbff_path = lt.get_path() + '/data/genomes/genomes_ncbi/ATCC13985/GCF_006322025_2_ASM632202v2_genomic.gbff'
        elif taxon == 'KBS0715':
            df_gbff_path = lt.get_path() + '/data/genomes/genomes_ncbi/KBS0715/GCF_005938005_2_ASM593800v2_genomic.gbff'
        else:
            continue
        for record in SeqIO.parse(df_gbff_path, "genbank"):
            for feature in record.features:
                if feature.type == 'CDS':
                    if feature.qualifiers['locus_tag'][0] in gene_names:
                        location = str(feature.location)
                        dist_dict[feature.qualifiers['locus_tag'][0]] = int(re.split(r"[:[']+", location)[1])

        dist_df = pd.DataFrame(list(dist_dict.items()), columns=['Name', 'Position'])
        dist_df.index = dist_df.Name
        del dist_df['Name']
        df_merge = df_mult.merge(dist_df, left_index=True, right_index=True)
        ax = fig.add_subplot(2, 2, i+1)
        ax.scatter(df_merge.Position.values, df_merge.Multiplicity.values, c='#175ac6', marker = 'o', s = 70, \
            edgecolors='#244162', linewidth = 0.6, alpha = 0.5, zorder=2)#, edgecolors='none')
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.1e'))

        ax.title.set_text(latex_dict[taxon])

    fig.text(0.5, 0.02, 'Genome position', ha='center', fontsize=16)
    fig.text(0.02, 0.5, 'Gene multiplicity, ' + '$m$', va='center', rotation='vertical', fontsize=16)

    fig.savefig(lt.get_path() + '/figs/dist_mult.png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()




def plot_birth_per_death_vs_shape():
    df_gen = pd.read_csv(lt.get_path() + '/data/breseq/genetic_diversity.txt', sep = '\t')
    df_gen = df_gen.rename(columns={'Species': 'strain'})
    df_weib = pd.read_csv(lt.get_path() + '/data/demography/weibull_results_clean.csv', sep = ',')
    df_merged = df_weib.merge(df_gen, on=['strain','rep'])
    taxa = list(set(df_merged.strain.to_list()))

    fig = plt.figure()

    plt.scatter(df_merged.max_N_mut.values, df_merged.alpha.values, color='k', linestyle='--', linewidth=2, zorder=2)

    plt.xlim(min(df_merged.max_N_mut.values), max(df_merged.max_N_mut.values))
    plt.ylim(min(df_merged.alpha.values), 1)

    plt.xscale('log',basex=10)
    plt.yscale('log',basey=10)
    plt.xlabel('births-per-deaths (hrs.)', fontsize = 12)
    plt.ylabel('Shape paramter, ' r'$k$', fontsize = 12)
    fig.savefig(lt.get_path() + '/figs/birth_per_death_vs_shape.png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()




def plot_ode():
    ts = np.linspace(0, 1000, 10000)
    N_0 = int(1e9)
    P0 = [N_0, 0, 0]
    colors = ['#FF6347', '#FFA500', '#87CEEB']
    labels = [r'$N(0)\cdot d = 10^{6}$', r'$N(0)\cdot d = 10^{7}$', r'$N(0)\cdot d = 10^{8}$']
    fig = plt.figure()

    for i, d in enumerate([ 0.001, 0.01, 0.1 ]):
        Ps = odeint(dP_dt, P0, ts, args=(d,))
        # Ps = odeint(dP_dt, P0, ts)
        N = Ps[:,0]

        plt.plot(ts, N, "-", c = colors[i], label=labels[i])

    plt.legend(loc='lower left', prop={'size': 8})

    plt.yscale('log',basey=10)

    plt.xlabel('Time, ' + r'$t$' , fontsize = 14)
    plt.ylabel('Number of cells, ' + r'$N(t)$' , fontsize = 14)

    fig.savefig(lt.get_path() + '/figs/ode_example.pdf', format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
