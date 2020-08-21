from __future__ import division
import os
import ltde_tools as lt

to_ignore = ['strain', 'avg']

spo2e_time = None
wt_time = None

spo2e_dilution_V = None
spo2e_dilution_S = None

wt_dilution_V = None
wt_dilution_S = None

output_file = open('/Users/WRShoemaker/GitHub/LTDE/data/demography/spo0IIE_assay.csv',"w")

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
