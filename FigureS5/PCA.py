import os,sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from collections import Counter

synFile = 'Grass_syntenyList.csv'
spList = ['sorghum3', 'maize1_v4', 'maize2_v4', 'paspalum']
dfsynlist = pd.read_csv(synFile, usecols = spList)
dfsynlist = dfsynlist.melt(id_vars=("sorghum3","paspalum"),value_name='maize', var_name='subgenome')
dfsynlist = dfsynlist.loc[dfsynlist['paspalum']!='No Gene']
dfsynlist = dfsynlist.loc[dfsynlist['maize']!='No Gene']
dfsynlist

dfpv = pd.read_csv('Pv_tpm.csv',sep='\t',index_col = 0)
dfzm = pd.read_csv('Zm_tpm.csv',sep='\t',index_col = 0)
dfzm.index = [x.split('_')[0] for x in dfzm.index]
dfsb = pd.read_csv('Sb_tpm.csv',sep='\t',index_col = 0)
dfsb.index = [x.rsplit('.',1)[0] for x in dfsb.index]
dfpv = dfpv.loc[dfsynlist['paspalum']]
dfzm = dfzm.loc[dfsynlist['maize'].values]
dfsb = dfsb.loc[dfsynlist['sorghum3'].values]
dfpv['mean'] = dfpv.mean(axis = 1)
dfzm['mean'] = dfzm.mean(axis = 1)
dfsb['mean'] = dfsb.mean(axis = 1)
dfpv = dfpv.loc[dfpv['mean']>50]
dfzm = dfzm.loc[dfzm['mean']>50]
dfsb = dfsb.loc[dfsb['mean']>50]

dfpv = np.log(dfpv)
dfzm = np.log(dfzm)
dfsb = np.log(dfsb)

import numpy as np
from sklearn.decomposition import PCA

def DoPCA(pcDf, species):
    pca = PCA(n_components=2)s
    expValue = pcDf.values
    principalComponents = pca.fit_transform(expValue)
    principalDf = pd.DataFrame(data = principalComponents, columns = ['PC1', 'PC2'])
    print(species, pca.explained_variance_ratio_)
    return principalDf

sns.set_style("darkgrid")
fig, axes = plt.subplots(1,3,figsize=(9,3),tight_layout=True)
###Paspalum 
condDict = {'F':'Full','N':'-N','P':'-P'}

pcPV =  dfpv[list(dfpv)[:-1]].transpose()
pcZm =  dfzm[list(dfzm)[:-1]].transpose()
pcSb =  dfsb[list(dfsb)[:-1]].transpose()

dfpcPv = DoPCA(pcPV, 'Paspalum')
dfpcZm = DoPCA(pcZm, 'Maize')
dfpcSb = DoPCA(pcSb, 'Sorghum')

dfpcPv['condition'] = [condDict[x[0]] for x in pcPV.index]
dfpcZm['condition'] = [condDict[x[0]] for x in pcZm.index]
dfpcSb['condition'] = [condDict[x[0]] for x in pcSb.index]

dfpcPv.to_csv('Paspalum_RNAseq_Lib_PCA.csv')
dfpcZm.to_csv('Maize_RNAseq_Lib_PCA.csv')
dfpcSb.to_csv('Sorghum_RNAseq_Lib_PCA.csv')

sns.scatterplot(x='PC1',y='PC2',hue = 'condition', data=dfpcPv, s=50, ax = axes[0] )
sns.scatterplot(x='PC1',y='PC2',hue = 'condition', data=dfpcZm, s=50, ax = axes[1] )
sns.scatterplot(x='PC1',y='PC2',hue = 'condition', data=dfpcSb, s=50, ax = axes[2] )



axes[0].set_xlabel('PC1 (0.79)',fontsize = 10)
axes[0].set_ylabel('PC2 (0.20)',fontsize = 10)
axes[0].set_title('P. vaginatum',fontstyle = 'italic', fontsize = 10)

axes[1].set_xlabel('PC1 (0.86)',fontsize = 10)
axes[1].set_ylabel('PC2 (0.10)',fontsize = 10)
axes[1].set_title('Z. mays',fontstyle = 'italic', fontsize = 10)

axes[2].set_xlabel('PC1 (0.66)',fontsize = 10)
axes[2].set_ylabel('PC2 (0.31)',fontsize = 10)
axes[2].set_title('S. bicolor', fontstyle = 'italic', fontsize = 10)

plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.legend(fontsize=10)

plt.savefig('PCA.svg')
