from __future__ import division
import os
from biom import parse_table
from biom import load_table
import pandas as pd
import numpy as np
import  matplotlib.pyplot as plt
from scipy import stats

mydir = os.path.expanduser("~/GitHub/LTDE/")


def clean_biom(genus = 'Bacillus'):
    #out_directory = mydir + 'data/EMP/EMP-clean/processed_data/' + genus
    #if not os.path.exists(out_directory):
    #    os.makedirs(out_directory)
    df_list = []
    directory = mydir + 'data/EMP/EMP-public/processed_data'
    for filename in os.listdir(directory):
        if filename.endswith(".biom"):
            print(filename)
            biom_path = os.path.join(directory, filename)
            table = load_table(biom_path)
            # get options for biom
            # dir(table)
            tax = table.metadata_to_dataframe('observation')
            tax_genus = tax[tax['taxonomy_5'].str.contains("g__" + genus)]
            genus_otus = tax_genus.index.values
            genus_otus_ = [str(x) for x in genus_otus]
            genus_otus_str = "|".join(genus_otus_)

            df = (table.to_dataframe())
            N = df.sum(0)
            df_genus = df[df.index.str.contains(genus_otus_str)]
            N_genus = df_genus.sum(0)
            N_merge = pd.concat([N, N_genus], axis=1)
            N_merge.columns = ['N', 'N_genus']
            N_merge = N_merge[N_merge.N_genus != 0]
            df_list.append(N_merge)

        else:
            continue

    #OUT_file_path = out_directory + '/' + filename.split('.')[0] + '.txt'
    dfs_concat = pd.concat(df_list)
    OUT_file_path = mydir + 'data/EMP/EMP-clean/processed_data/' + genus + '.txt'
    dfs_concat.to_csv(OUT_file_path, sep = '\t', index = True)



def bacillus_fig():
    df_path = mydir + 'data/EMP/EMP-clean/processed_data/Bacillus.txt'
    df = pd.read_csv(df_path, sep = '\t')
    df = df[df.N_genus >= 10]
    x = np.log10(df.N.values)
    y = np.log10(df.N_genus.values)

    fig = plt.figure()

    plt.scatter(x, y, c='#87CEEB', marker='o', alpha = 0.4)
    #plt.plot([-5, 11],[-5, 11], '--', color = 'dimgrey')

    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    predict_y = intercept + slope * x
    pred_error = y - predict_y
    degrees_of_freedom = len(x) - 2
    residual_std_error = np.sqrt(np.sum(pred_error**2) / degrees_of_freedom)
    #plt.plot(x, predict_y, 'k-', lw = 1.5, c = 'black', label = 'm = ' + str(slope) )

    print(slope)
    print(r_value ** 2)

    plt.xlabel('Total site abundance, log10', fontsize = 14)
    plt.ylabel('Bacillus abundance, log10', fontsize = 14)

    fig.tight_layout()
    fig_name = mydir + 'figs/bacillus_N.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()


def get_spore_data():
    bac_path = mydir + 'data/bacillus_test/evo12597-sup-0004-Table-S2.txt'
    bac = pd.read_csv(bac_path, sep = '\t')
    bac = bac.drop(bac.index[197])
    def split_data(x):
        if x:
            s = x.split(' ')
            genus = s[0].strip()
            if s[1].strip() == 'sp.':
                species = s[2].strip()
            else:
                species = s[1].strip()
        else:
            genus = None
            species = None
        return pd.Series([genus, species], index=['genus','species'])

    bac[['genus','species']] = bac['Taxon Name'].apply(split_data)
    spore_y = bac.loc[bac['Spore-Forming?'] == 'Y']
    spore_n = bac.loc[bac['Spore-Forming?'] == 'N']
    spore_y_zip =  list(zip(spore_y.genus, spore_y.species))
    spore_n_zip =  list(zip(spore_n.genus, spore_n.species))

    df_list = []
    directory = mydir + 'data/EMP/EMP-public/processed_data'
    for filename in os.listdir(directory):
        if filename.endswith(".biom"):
            print(filename)
            biom_path = os.path.join(directory, filename)
            table = load_table(biom_path)
            tax = table.metadata_to_dataframe('observation')

            tax_sporeY_list = []
            for i in spore_y_zip:
                spore_tax_i = tax[(tax["taxonomy_5"].str.contains("g__" + i[0]) & (tax["taxonomy_6"].str.contains("s__" + i[1])))]
                tax_sporeY_list.append(spore_tax_i)
            tax_sporeN_list = []
            for j in spore_n_zip:
                spore_tax_j = tax[(tax["taxonomy_5"].str.contains("g__" + j[0]) & (tax["taxonomy_6"].str.contains("s__" + j[1])))]
                tax_sporeN_list.append(spore_tax_j)

            tax_sporeN = pd.concat(tax_sporeN_list)
            tax_sporeN_otus = tax_sporeN.index.values
            tax_sporeN_otus_ = [str(x) for x in tax_sporeN_otus]
            tax_sporeN_otus_str = "|".join(tax_sporeN_otus_)

            tax_sporeY = pd.concat(tax_sporeY_list)
            tax_sporeY_otus = tax_sporeY.index.values
            tax_sporeY_otus_ = [str(x) for x in tax_sporeY_otus]
            tax_sporeY_otus_str = "|".join(tax_sporeY_otus_)

            df = table.to_dataframe()
            N = df.sum(0)
            S = df.astype(bool).sum(axis=0)

            df_sporeN = df[df.index.str.contains(tax_sporeN_otus_str)]
            N_sporeN = df_sporeN.sum(0)
            S_sporeN = df_sporeN.astype(bool).sum(axis=0)
            df_sporeY = df[df.index.str.contains(tax_sporeY_otus_str)]
            N_sporeY = df_sporeY.sum(0)
            S_sporeY = df_sporeY.astype(bool).sum(axis=0)

            N_merge = pd.concat([N, S, N_sporeN, S_sporeN, N_sporeY, S_sporeY], axis=1)
            N_merge.columns = ['N','S','N_sporeN','S_sporeN', 'N_sporeY', 'S_sporeY']
            N_merge = N_merge[(N_merge.N_sporeN != 0) | (N_merge.N_sporeY != 0)]
            df_list.append(N_merge)

    dfs_concat = pd.concat(df_list)
    OUT_file_path = mydir + 'data/EMP/EMP-clean/processed_data/spore_counts.txt'
    dfs_concat.to_csv(OUT_file_path, sep = '\t', index = True)


def get_spore_metadata(env_variable = 'temp'):
    df_path = mydir + 'data/EMP/EMP-clean/processed_data/spore_counts.txt'
    df = pd.read_csv(df_path, sep = '\t', index_col = 0)
    sites = df.index.values
    sites_ = [str(x) for x in sites]
    sites_str = "|".join(sites_)
    directory = mydir + 'data/EMP/EMP-public/templates'
    temps_list = []
    for filename in os.listdir(directory):
        if filename.endswith(".txt")  & ("_prep_" not in filename):
            print(filename)
            temp_path = os.path.join(directory, filename)
            temp = pd.read_csv(temp_path, sep = '\t', index_col = 0)
            temp.index = temp.index.map(str)
            if env_variable == 'temp':
                if ('temp' in temp) or ('temperature_deg_c' in temp):
                    if ('temp' in temp) & ('temperature_deg_c' in temp):
                        temp = temp.loc[:,temp.columns.isin(["temp", "temperature_deg_c"])]
                        temp[['temp', 'temperature_deg_c']] = temp[['temp', 'temperature_deg_c']].apply(pd.to_numeric, errors='coerce')
                    else:
                        temp = temp.loc[:,temp.columns.isin(["temp"])]
                        temp[['temp']] = temp[['temp']].apply(pd.to_numeric, errors='coerce')
                    temp = temp[np.isfinite(temp['temp'])]
                    temp = temp[pd.notnull(temp['temp'])]
                    temp = temp[['temp']]
                    #temp = temp[['temp']]
                    #print(temp)
                else:
                    continue

            elif env_variable == 'ph':
                if 'ph' in temp:
                    temp = temp.loc[:,temp.columns.isin(["ph"])]
                    temp[['ph']] = temp[['ph']].apply(pd.to_numeric, errors='coerce')
                    temp = temp[np.isfinite(temp['ph'])]
                    temp = temp[pd.notnull(temp['ph'])]
                    temp = temp[['ph']]
                    if temp.shape[0] == 0:
                        continue
                else:
                    continue
            else:
                "Choose an environmental variable"
                break

            temp_keep = temp[temp.index.str.contains(sites_str)]
            if temp_keep.shape[0] == 0:
                continue
            temps_list.append(temp_keep)
    temps_concat = pd.concat(temps_list)
    #print(temps_concat)
    # merge the environmental df and diversity df
    df_temp = pd.merge(df, temps_concat, left_index=True, right_index=True)
    df_temp_path = mydir + 'data/EMP/EMP-clean/processed_data/spore_counts_' + env_variable + '.txt'
    df_temp.to_csv(df_temp_path, sep = '\t', index = True)


def make_spore_plots_temp():
    df_path = mydir + 'data/EMP/EMP-clean/processed_data/spore_counts_temp.txt'
    df = pd.read_csv(df_path, sep = '\t', index_col = 0)
    df_sporeY = df[df.S_sporeY != float(0)]
    #df_sporeY = df_sporeY.loc[df_sporeY['S_sporeY'] < 500]

    df_sporeN = df[df.S_sporeN != float(0)]
    #df_sporeN = df_sporeN.loc[df_sporeN['S_sporeN'] < 500]
    # make plots twice
    x_Y = df_sporeY.temp.values
    y_Y = np.log10(df_sporeY.S_sporeY.values)

    fig = plt.figure()
    plt.scatter(x_Y, y_Y, c='#87CEEB', marker='o', alpha = 0.4)

    plt.xlabel('Temp', fontsize = 14)
    plt.ylabel('Richness, spore formners', fontsize = 14)
    plt.xlim(-20,100)
    plt.ylim(0,3.7)

    fig.tight_layout()
    fig_name = mydir + 'figs/spore_temp.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()


    x_N = df_sporeN.temp.values
    y_N = np.log10(df_sporeN.S_sporeN.values)

    fig = plt.figure()
    plt.scatter(x_N, y_N, c='#87CEEB', marker='o', alpha = 0.4)

    plt.xlabel('Temp', fontsize = 14)
    plt.ylabel('Richness, non-spore formners', fontsize = 14)
    plt.xlim(-20,100)
    plt.ylim(0,3.7)

    fig.tight_layout()
    fig_name = mydir + 'figs/non_spore_temp.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()




def make_spore_plots_ph():
    df_path = mydir + 'data/EMP/EMP-clean/processed_data/spore_counts_ph.txt'
    df = pd.read_csv(df_path, sep = '\t', index_col = 0)
    df_sporeY = df[df.S_sporeY != float(0)]
    #df_sporeY = df_sporeY.loc[df_sporeY['S_sporeY'] < 500]

    df_sporeN = df[df.S_sporeN != float(0)]
    #df_sporeN = df_sporeN.loc[df_sporeN['S_sporeN'] < 500]
    # make plots twice
    x_Y = df_sporeY.ph.values
    y_Y = np.log10(df_sporeY.S_sporeY.values)

    fig = plt.figure()
    plt.scatter(x_Y, y_Y, c='#87CEEB', marker='o', alpha = 0.4)

    plt.xlabel('ph', fontsize = 14)
    plt.ylabel('Richness, spore formners, log10', fontsize = 14)
    plt.xlim(2,12)
    plt.ylim(0,3.7)

    fig.tight_layout()
    fig_name = mydir + 'figs/spore_ph.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()


    x_N = df_sporeN.ph.values
    y_N = np.log10(df_sporeN.S_sporeN.values)
    testtt = df_sporeN.loc[df_sporeN['S_sporeN'] > 500]
    print(testtt)

    fig = plt.figure()
    plt.scatter(x_N, y_N, c='#87CEEB', marker='o', alpha = 0.4)

    plt.xlabel('ph', fontsize = 14)
    plt.ylabel('Richness, non-spore formers, log10', fontsize = 14)
    plt.xlim(2,12)
    plt.ylim(0,3.7)

    fig.tight_layout()
    fig_name = mydir + 'figs/non_spore_ph.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()


# get sites
# merge sites from all templates into single df
# select sites that contain
# pH, temperature, salinity, and oxygen
# latitude and elevation

#clean_biom()
#bacillus_fig()
#get_spore_data()
#get_spore_metadata()
#get_spore_metadata(env_variable = 'temp')
#get_spore_metadata(env_variable = 'ph')
make_spore_plots_temp()
#make_spore_plots_ph()
