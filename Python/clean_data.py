from __future__ import division
import ltde_tools as lt
import glob, re, os, subprocess, math
import pandas as pd

def get_iRep():
    directory = os.fsencode(lt.get_path() + '/data/bwa_sam')
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if (filename.endswith('C1.sam') == True) or (filename.endswith('C2.sam') == True):
            continue
        if filename.endswith('.sam'):
            strain = re.split(r'[.-]+', filename)[0]
            print(filename)
            fna_path = glob.glob(lt.get_path() + '/data/genomes/*/' + strain + '/G-Chr1.fna')[0]
            out_file = lt.get_path() + '/data/iRep/' + filename.split('.')[0]
            sam = os.path.join(str(directory, 'utf-8'), filename)
            subprocess.call(['iRep', '-f', fna_path, '-s', sam, '-o', str(out_file)])


def clean_iRep():
    to_remove = ['KBS0705', 'KBS0706']
    directory = os.fsencode(lt.get_path() + '/data/iRep')
    df_out = open(lt.get_path() + '/data/iRep_clean.txt', 'w')
    header = ['Sample', 'strain', 'rep' ,'iRep']
    df_out.write('\t'.join(header) + '\n')
    iRep_corrected_dict = {}
    iRep_uncorrected_dict = {}
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith('.tsv'):
            iRep_path = sam = os.path.join(str(directory, 'utf-8'), filename)
            strain = re.split(r'[.-]+', filename)[0]
            strain_rep = re.split(r'[.]+', filename)[0]
            if strain in to_remove:
                continue
            if 'WA' in strain_rep:
                strain_rep = 'KBS0711-5'
            elif 'WB' in strain_rep:
                strain_rep = 'KBS0711-6'
            elif 'WC' in strain_rep:
                strain_rep = 'KBS0711-7'
            else:
                strain_rep = strain_rep[:-1] + str(lt.rename_rep()[strain_rep[-1]])
            for i, line in enumerate(open(iRep_path, 'r')):
                if i == 2:
                    last_item = line.strip().split()[-1]
                    if last_item == 'n/a':
                        iRep_corrected = float('nan')
                    else:
                        iRep_corrected = float(last_item)
                    iRep_corrected_dict[strain_rep] = [iRep_corrected]
                elif i == 6:
                    iRep_uncorrected = float(line.strip().split()[-1])
                    iRep_uncorrected_dict[strain_rep] = [iRep_uncorrected]
    for key, value in iRep_corrected_dict.items():
        value.extend(iRep_uncorrected_dict[key])
    for key, value in iRep_corrected_dict.items():
        if value[1] > 11:
            continue

        if math.isnan(value[0]) == True:
            iRep = value[1]
        else:
            iRep = value[0]
        out_line = [key, key.split('-')[0], key.split('-')[1], str(iRep)]
        df_out.write('\t'.join(out_line) + '\n')


    df_out.close()


def clean_COGs():
    directory = os.fsencode(lt.get_path() + '/data/COGs')
    cog_dict = {}
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith('_hmms.hits.txt'):
            cog_path = sam = os.path.join(str(directory, 'utf-8'), filename)
            df = pd.read_csv(cog_path, sep = ',')
            strain = filename.split('_')[0]
            cog_dict[strain] = {}
            cog_list = [x for x in df.gene_hmm_id.tolist() if x != '-']
            for cog in cog_list:
                cog_dict[strain][cog] = 1

    df_cogs = pd.DataFrame.from_dict(cog_dict)
    df_cogs = df_cogs.fillna(0)
    #df_cogs = df_cogs[(df_cogs.T != 1).any()]
    #df_cogs = df_cogs[(df_cogs != 1).any()]
    #print(df_cogs)
    df_cogs = df_cogs[(df_cogs.T != 1).any()].T
    df_out = lt.get_path() + '/data/COGs/cog_by_genome.txt'
    df_cogs.to_csv(df_out, sep = '\t', index = True)



#clean_iRep()
clean_COGs()
