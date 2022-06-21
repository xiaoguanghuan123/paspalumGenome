import pandas as pd
import os 
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rcParams

dfsynList = pd.read_csv('SbZmPvSynListUsedForDEG.csv',index_col = 'CommonID')
dfPF = pd.DataFrame(index = dfsynList.index)


dfsbPF = pd.read_csv('Sorghum_PF_adjp0.05.csv',index_col = 0)
dfzmPF = pd.read_csv('Maize_PF_adjp0.05.csv',index_col = 0)
dfpvPF = pd.read_csv('Paspalum_PF_adjp0.05.csv',index_col = 0)

dfpvPF.loc[dfpvPF.log2FoldChange >= 5, 'log2FoldChange' ] = 5
dfpvPF.loc[dfpvPF.log2FoldChange <= -5, 'log2FoldChange'] = -5
dfzmPF.loc[dfzmPF.log2FoldChange >= 5, 'log2FoldChange' ]= 5
dfzmPF.loc[dfzmPF.log2FoldChange <= -5, 'log2FoldChange' ]= -5
dfsbPF.loc[dfsbPF.log2FoldChange >= 5, 'log2FoldChange' ]= 5
dfsbPF.loc[dfsbPF.log2FoldChange <= -5, 'log2FoldChange' ]= -5

dfpvPF = dfpvPF.loc[(abs(dfpvPF.log2FoldChange) >= 1)]
dfzmPF = dfzmPF.loc[(abs(dfzmPF.log2FoldChange) >= 1)]
dfsbPF = dfsbPF.loc[(abs(dfsbPF.log2FoldChange) >= 1)]

def lg2fc(df,agene,clist):
    if agene in df.index:
        lg2 = df.loc[agene].log2FoldChange
        clist.append(lg2)
    else: clist.append(0)

sbList, zmList, pvList = [], [], []
for agene in dfPF.index:
    sbGene = dfsynList.loc[agene]['sorghum3']
    zmGene = dfsynList.loc[agene]['maize4']
    pvGene = dfsynList.loc[agene]['paspalum3']+'.g'
    lg2fc(dfsbPF, sbGene, sbList)
    lg2fc(dfzmPF, zmGene, zmList)
    lg2fc(dfpvPF, pvGene, pvList)
dfPF['paspalum'] = pvList
dfPF['maize'] = zmList
dfPF['sorghum'] = sbList

for agene in dfPF.index:
    vCount = dfPF.loc[agene].value_counts()
    if 0 in vCount and vCount[0] == 3: dfPF = dfPF.drop(agene, axis = 0)
        
print(dfPF.info())
rcParams['font.family'] = 'sans-serif'
rcParams['font.size'] = 10
fig, ax = plt.subplots(figsize=(1,10),tight_layout = True)
sns.heatmap(dfPF, annot = False, cmap = 'coolwarm', cbar = False, center=0, ax=ax) #cbar_kws={"orientation": "horizontal"}
ax.set_yticks([])
ax.set_yticklabels([])
ax.set_ylabel('')
ax.set_xticklabels(['Pv','Zm','Sb'],rotation = 360, fontstyle = 'italic')
ax.set_title('-P/Full')

plt.savefig('FP_Heatmap_1.svg')
plt.savefig('FP_Heatmap_1.png',dpi=500)
