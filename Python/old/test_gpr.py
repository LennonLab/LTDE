import os
import numpy as np
import pandas as pd
from fitderiv import fitderiv
import matplotlib.pyplot as plt

df = pd.read_csv(os.path.expanduser("~/GitHub/LTDE") + '/data/demography/longtermdormancy_20190528_nocomments.csv', sep = ',')
# KBS0721 rep
df['N'] = (df['Colonies']+1) * (1000 / df['Inoculum']) * (10 ** (df['Dilution'] ))

df['Dormstart_date'] =  pd.to_datetime(df['Dormstart_date'], format='%d-%b-%y')
df['Firstread_date'] =  pd.to_datetime(df['Firstread_date'], format='%d-%b-%y')
df['Days'] = df['Firstread_date'].sub(df['Dormstart_date'], axis=0)
df['Days'] = df['Days'].dt.days.astype('int')

df_taxon = df[(df["Strain"] == 'KBS0706')]
#df_taxon = df_taxon.sort_values(by ='Days')
days_list = np.sort(list(set(df_taxon.Days.tolist())))
print(df_taxon)
#arr = np.empty((0,5), float)
lists = []
for days in days_list:
    df_taxon_days =  df_taxon[(df_taxon["Days"] == days)]
    df_taxon_days = df_taxon_days.sort_values(by ='Rep')
    print([days] + df_taxon_days.N.tolist())
    if len([days] + df_taxon_days.N.tolist()) != 5:
        continue
    #data_list.append(np.asarray([days] + df_taxon_days.N.tolist()))
    lists.append([days] + df_taxon_days.N.tolist())
    #data_list.append()

arrayyy = np.asarray(lists)
print(arrayyy)

t= arrayyy[:, 0]
od= arrayyy[:, 1:]

# call fitting algorithm
# try help(fitderiv) for more information
q = fitderiv(t, od, cvfn= 'sqexp', stats= True, esterrs= True)

print(q.ddf)


plt.figure()
plt.subplot(3,1,1)
q.plotfit('f')
plt.subplot(3,1,2)
q.plotfit('df')
plt.subplot(3,1,3)
q.plotfit('ddf')

plt.show()

#print(q.ds['max df'])
