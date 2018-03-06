from __future__ import division
import os, re, shutil
from subprocess import call


# open archive-GSF1046-22Apr2016.tar.gz into reads_raw in bash
# nohup tar -zxvf archive-GSF1046-22Apr2016.tar.gz -C /N/dc2/projects/muri2/Task2/LTDE/data/reads_raw


def get_index_dict():
    index_path = '/N/dc2/projects/muri2/Task2/LTDE/bin/GSF1046_libraries.txt'
    index_dict = {}
    for line in open(index_path, "r"):
        line_clean = line.strip().split('\t')
        if line_clean[0] == 'Sample':
            continue
        index_dict[line_clean[0]] = line_clean[1]
    return(index_dict)


def rename():
    index_dict = get_index_dict()
    lib_path = '/N/dc2/projects/muri2/Task2/LTDE/data/reads_raw'
    for lib in os.listdir(lib_path):
        #lib_rename =
        lib_sample = lib.split('_')[0]
        indexes = index_dict[lib_sample]
        lib_rename_split = lib.split('.', 1)
        lib_rename = lib_rename_split[0] + '_' + indexes + '.' + lib_rename_split[1]
        out = '/N/dc2/projects/muri2/Task2/LTDE/data/reads_raw_rename/' + lib_rename
        shutil.copy2(lib_path + '/' + lib, out)


rename()
#    for i, x in enumerate(list(os.walk('/N/dc2/projects/muri2/Task2/PoolPopSeq/data/' + run))):
#if i  == 0:
#    continue
#path = x[0]
#for y in x[2]:
#    if y == 'SampleSheet.csv':
#        continue
#    if run == '161213':
#        y_split = y.split('_', 1)
#        y_out = y_split[0] + '-100_' + y_split[1].split('.', 1)[0] + '_' + \
#            run + '.' + y_split[1].split('.', 1)[1]
#    elif run == '170303':
#        y_out = y.split('.', 1)[0] + '_' + run + '.' + y.split('.', 1)[1]
#    elif run == '170721':
#        y_rep = y.replace('-D10', '')
#        y_rep2 = y_rep.replace('OF3', '0F3')
#        y_out = 'L' + y_rep2.split('.', 1)[0] + '_' + run + '.' + y_rep2.split('.', 1)[1]
#    elif run == '170623':
#        y_split = y.split('_', 4)
#        y_out = 'L' + y_split[1] + y_split[2] + '-' + y_split[3][1:] \
#            + '_' + y_split[4].split('.', 1)[0] + '_' + run + '.' + y_split[4].split('.', 1)[1]
#    y_out_path1 = '/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_raw/' \
#            + 'D' +  re.split('-|_',y_out)[1] + '/'
#    if os.path.exists(y_out_path1) != True:
#        os.makedirs(y_out_path1)
#    y_out_path2 = y_out_path1 + 'Sample_' + y_out.split('_', 1)[0] + '/'
#    if os.path.exists(y_out_path2) != True:
#        os.makedirs(y_out_path2)
#    y_out_path = y_out_path2 + y_out
#    y_path = x[0] + '/' + y
#    call(["cp", y_path, y_out_path])
#    print y_path
#    print y_out_path
