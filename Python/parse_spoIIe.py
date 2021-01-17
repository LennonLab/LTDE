from __future__ import division
import ltde_tools as lt
import glob, re, os, subprocess, math, json, copy
import pandas as pd
import numpy as np
from collections import Counter
from Bio import SeqIO
import _pickle as pickle
from decimal import Decimal
import statsmodels.stats.multitest as mt
from scipy.stats import t

import matplotlib.pyplot as plt


def clean_SpoIIE_GC_data_new():

    to_ignore = ['strain', 'avg']

    spo2e_time = None
    wt_time = None

    spo2e_dilution_V = None
    spo2e_dilution_S = None

    wt_dilution_V = None
    wt_dilution_S = None

    output_file = open(lt.get_path() + '/data/demography/spo0IIE_assay.csv',"w")

    spo2e_header = ",".join(['strain', 'replicate', 'NT', 'HT', 'hours', 'dilution_V', 'dilution_S' ])

    output_file.write(spo2e_header)
    output_file.write('\n')

    for line_idx, line in enumerate(open('/Users/WRShoemaker/Desktop/SpoIIE_GC.csv')):
        line = line.strip().split(',')
        if (line_idx in list(range(4))) or (line[0] in to_ignore):
            continue

        spo2e_rep = int(line[0].split(' ')[1])
        wt_rep = int(line[8].split(' ')[1])

        if len(line[3]) > 0:
            spo2e_time = float(line[3].replace("hr",""))

        if len(line[11]) > 0:
            wt_time = float(line[3].replace("hr",""))



        if (len(line[1]) > 0):
            spo2e_NT = int(line[1])
            if len(line[2]) == 0:
                spo2e_HT = 0
            else:
                spo2e_HT = int(line[2])

            if (len(line[5]) > 0) and "CFU's" not in line[5]:
                spo2e_dilution_V = abs(int(line[5].split('^')[1]))
                spo2e_dilution_S = abs(int(line[6].split('^')[1]))


            spo2e_line = ",".join(['SpoIIE', str(spo2e_rep), str(spo2e_NT), str(spo2e_HT),  str(spo2e_time), str(spo2e_dilution_V), str(spo2e_dilution_S) ])

            output_file.write(spo2e_line)
            output_file.write('\n')

        if (len(line[9]) > 0):
            wt_NT = int(line[9])
            if len(line[10]) == 0:
                wt_HT = 0
            else:
                wt_HT = int(line[10])

            if (len(line[13]) > 0) and "CFU's" not in line[13]:
                wt_dilution_V = abs(int(line[13].split('^')[1]))
                wt_dilution_S = abs(int(line[14].split('^')[1]))

                if ( wt_time > 17*24) :
                    wt_dilution_V = 4
                    wt_dilution_S = 4

            wt_line = ",".join(['wt', str(wt_rep), str(wt_NT), str(wt_HT),  str(wt_time), str(wt_dilution_V), str(wt_dilution_S) ])

            output_file.write(wt_line)
            output_file.write('\n')



    output_file.close()






def clean_SpoIIE_GC_data():
    file = open(lt.get_path() + '/data/demography/SpoIIE_GC_assay.csv', "w")
    header = ','.join(['strain', 'replicate', 'NT', 'HT', 'hours', 'dilution_V', 'dilution_S']) + '\n'
    file.write(header)

    to_ignore = ['Bacillus', 'strains', 'medium', 'temp', '', 'strain', 'avg', 'temp - 30C']

    time_spo_all = None
    time_wt_all = None

    dilution_V_spo_all = None
    dilution_V_wt_all = None
    dilution_S_spo_all = None
    dilution_S_wt_all = None

    for line in open(lt.get_path() + '/data/demography/SpoIIE_GC.csv', 'r'):

        line_split = line.strip().split(',')

        if ('strains' in line_split[0]) or (line_split[1] == '') or ('SpoIIE' not in line_split[0]):
            continue

        print(line_split)


        #print(line_split[8])
        type_spo = 'SpoIIE'
        rep_spo = line_split[0].split(' ')[1]
        type_wt = 'wt'
        rep_wt = line_split[8].split(' ')[1]

        NT_spo = line_split[1]
        NT_wt = line_split[9]

        HT_spo = line_split[2]
        HT_wt = line_split[10]



        time_spo = line_split[3].replace('hr', '')
        time_wt = line_split[11].replace('hr', '')

        if time_spo != '':
            time_spo_all = time_spo
            time_wt_all = time_wt

        dilution_V_spo = line_split[5].replace('10^', '')
        dilution_V_spo = dilution_V_spo.replace('-', '')
        if (dilution_V_spo != '') and (dilution_V_spo != "CFU's"):
            dilution_V_spo_all = dilution_V_spo

        dilution_V_wt = line_split[13].replace('10^', '')
        dilution_V_wt = dilution_V_wt.replace('-', '')
        if (dilution_V_wt != '') and (dilution_V_wt != "CFU's"):
            dilution_V_wt_all = dilution_V_wt

        dilution_S_spo = line_split[6].replace('10^', '')
        dilution_S_spo = dilution_S_spo.replace('-', '')
        if (dilution_S_spo != '') and (dilution_S_spo != "CFU's"):
            dilution_S_spo_all = dilution_S_spo

        dilution_S_wt = line_split[14].replace('10^', '')
        dilution_S_wt = dilution_S_wt.replace('-', '')
        if (dilution_S_wt != '') and (dilution_S_wt != "CFU's"):
            dilution_S_wt_all = dilution_S_wt

        line_spo = ','.join([type_spo, rep_spo, NT_spo, HT_spo, time_spo_all, dilution_V_spo_all, dilution_S_spo_all ])
        line_wt = ','.join([type_wt, rep_wt, NT_wt, HT_wt, time_wt_all, dilution_V_wt_all, dilution_S_wt_all ])

        file.write(line_spo+'\n')
        file.write(line_wt+'\n')
        # one bad data point, don't know why don't care
        #if '4.30E+06' not in line_wt:
        #    file.write(line_wt+'\n')

    file.close()



def clean_SpoIIE_DC_data():
    file = open(lt.get_path() + '/data/demography/SpoIIE_DC_assay.csv', "w")
    header = ','.join(['strain', 'replicate', 'NT', 'HT', 'hours', 'dilution_V', 'dilution_S']) + '\n'
    file.write(header)

    to_ignore = ['Bacillus', 'strains', 'medium', 'temp', '', 'strain', 'avg', 'temp - 30C']

    time_spo_all = None
    time_wt_all = None

    dilution_V_spo_all = None
    dilution_V_wt_all = None
    dilution_S_spo_all = None
    dilution_S_wt_all = None

    for line in open(lt.get_path() + '/data/demography/SpoIIE_DC.csv', 'r'):

        line_split = line.strip().split(',')

        if ('strains' in line_split[0]) or (line_split[1] == '') or ('SpoIIE' not in line_split[0]):
            continue


        #print(line_split[8])
        type_spo = 'SpoIIE'
        rep_spo = line_split[0].split(' ')[1]
        type_wt = 'wt'
        rep_wt = line_split[8].split(' ')[1]

        NT_spo = line_split[1]
        NT_wt = line_split[9]

        HT_spo = line_split[2]
        HT_wt = line_split[10]



        time_spo = line_split[3].replace('hr', '')
        time_wt = line_split[11].replace('hr', '')

        if time_spo != '':
            time_spo_all = time_spo
            time_wt_all = time_wt

        dilution_V_spo = line_split[5].replace('10^', '')
        dilution_V_spo = dilution_V_spo.replace('-', '')
        if (dilution_V_spo != '') and (dilution_V_spo != "CFU's"):
            dilution_V_spo_all = dilution_V_spo

        dilution_V_wt = line_split[13].replace('10^', '')
        dilution_V_wt = dilution_V_wt.replace('-', '')
        if (dilution_V_wt != '') and (dilution_V_wt != "CFU's"):
            dilution_V_wt_all = dilution_V_wt

        dilution_S_spo = line_split[6].replace('10^', '')
        dilution_S_spo = dilution_S_spo.replace('-', '')
        if (dilution_S_spo != '') and (dilution_S_spo != "CFU's"):
            dilution_S_spo_all = dilution_S_spo

        dilution_S_wt = line_split[14].replace('10^', '')
        dilution_S_wt = dilution_S_wt.replace('-', '')
        if (dilution_S_wt != '') and (dilution_S_wt != "CFU's"):
            dilution_S_wt_all = dilution_S_wt

        line_spo = ','.join([type_spo, rep_spo, NT_spo, HT_spo, time_spo_all, dilution_V_spo_all, dilution_S_spo_all ])
        line_wt = ','.join([type_wt, rep_wt, NT_wt, HT_wt, time_wt_all, dilution_V_wt_all, dilution_S_wt_all ])

        file.write(line_spo+'\n')
        # one bad data point, don't know why don't care
        if '4.30E+06' not in line_wt:
            file.write(line_wt+'\n')

    file.close()
