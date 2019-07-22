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
