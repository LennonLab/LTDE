import os
import pandas as pd
from sklearn.linear_model import ElasticNet
from sklearn.metrics import r2_score

mydir = os.path.expanduser("~/GitHub/LTDE/")

df = pd.read_csv(mydir + 'data/module_by_taxon.txt', sep = '\t', index_col = False)
df = df.T
df.columns = df.iloc[0]
df = df.drop(df.index[0])
df = df.reset_index(drop=False)
df.rename(columns={'index':'strain'}, inplace=True)

nunique = df.apply(pd.Series.nunique)
cols_to_drop = nunique[nunique == 1].index
df = df.drop(cols_to_drop, axis=1)

df_weibull = pd.read_csv(mydir + 'data/weibull_log_results.csv', sep = ',', index_col = False)

df_weibull.loc[df_weibull.strain == 'KBS0711W', 'strain'] = 'KBS0711'


df_merged = pd.merge(df, df_weibull, on='strain')

X = df_merged[df_merged.columns[df_merged.columns.to_series().str.contains('M0')]]
y = df_merged['Half_life']

X = X.values
y = y.values

n_samples = X.shape[0]
X_train, y_train = X[:n_samples // 2], y[:n_samples // 2]
X_test, y_test = X[n_samples // 2:], y[n_samples // 2:]


alpha = 0.1
enet = ElasticNet(alpha=alpha, l1_ratio=0.7)
y_pred_enet = enet.fit(X_train, y_train).predict(X_test)

r2_score_enet = r2_score(y_test, y_pred_enet)

#print("r^2 on test data : %f" % r2_score_enet)
