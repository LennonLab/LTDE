from __future__ import division
import ltde_tools as lt
import glob, re, os, subprocess

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
    directory = os.fsencode(lt.get_path() + '/data/iRep')
    df_out = open(lt.get_path() + '/data/iRep_clean.txt', 'w')
    header = ['Sample', 'Strain', 'Replicate' ,'iRep']
    df_out.write('\t'.join(header) + '\n')
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith('-100.txt'):
            bPTR_path = sam = os.path.join(str(directory, 'utf-8'), filename)
            for i, line in enumerate(open(bPTR_path, 'r')):
                if i == 0:
                    continue
                f_clean = filename.split('.')[0]
                f_clean_split = re.split(r'[-_]+', f_clean)
                out_line = [f_clean, f_clean_split[1][2],  f_clean_split[1][1],
                            f_clean_split[1][3], f_clean_split[2],  line.split()[-1]]
                df_out.write('\t'.join(out_line) + '\n')
                print(f_clean)
    df_out.close()


#get_iRep()
