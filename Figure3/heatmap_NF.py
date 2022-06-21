import pandas as pd
import os 
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rcParams

dfsynList = pd.read_csv('SbZmPvSynListUsedForDEG.csv',index_col = 'CommonID')
dfNF = pd.DataFrame(index = dfsynList.index)


dfsbNF = pd.read_csv('Sorghum_NF_adjp0.05.csv',index_col = 0)
dfzmNF = pd.read_csv('Maize_NF_adjp0.05.csv',index_col = 0)
dfpvNF = pd.read_csv('Paspalum_NF_adjp0.05.csv',index_col = 0)

dfpvNF.loc[dfpvNF.log2FoldChange >= 5, 'log2FoldChange' ] = 5
dfpvNF.loc[dfpvNF.log2FoldChange <= -5, 'log2FoldChange'] = -5
dfzmNF.loc[dfzmNF.log2FoldChange >= 5, 'log2FoldChange' ]= 5
dfzmNF.loc[dfzmNF.log2FoldChange <= -5, 'log2FoldChange' ]= -5
dfsbNF.loc[dfsbNF.log2FoldChange >= 5, 'log2FoldChange' ]= 5
dfsbNF.loc[dfsbNF.log2FoldChange <= -5, 'log2FoldChange' ]= -5


#dfpvNF = dfpvNF.loc[(abs(dfpvNF.log2FoldChange) >= 1) & (abs(dfpvNF.log2FoldChange) <= 5)]
#dfzmNF = dfzmNF.loc[(abs(dfzmNF.log2FoldChange) >= 1) & (abs(dfzmNF.log2FoldChange) <= 5)]
#dfsbNF = dfsbNF.loc[(abs(dfsbNF.log2FoldChange) >= 1) & (abs(dfsbNF.log2FoldChange) <= 5)]

def lg2fc(df,agene,clist):
    if agene in df.index:
        lg2 = df.loc[agene].log2FoldChange
        clist.append(lg2)
    else: clist.append(0)

sbList, zmList, pvList = [], [], []
for agene in dfNF.index:
    sbGene = dfsynList.loc[agene]['sorghum3']
    zmGene = dfsynList.loc[agene]['maize4']
    pvGene = dfsynList.loc[agene]['paspalum3']+'.g'
    lg2fc(dfsbNF, sbGene, sbList)
    lg2fc(dfzmNF, zmGene, zmList)
    lg2fc(dfpvNF, pvGene, pvList)
dfNF['paspalum'] = pvList
dfNF['maize'] = zmList
dfNF['sorghum'] = sbList

for agene in dfNF.index:
    vCount = dfNF.loc[agene].value_counts()
    if 0 in vCount and vCount[0] == 3: dfNF = dfNF.drop(agene, axis = 0)
        

rcParams['font.family'] = 'sans-serif'
rcParams['font.size'] = 10
fig, ax = plt.subplots(figsize=(1,3),tight_layout = True)
sns.heatmap(dfNF, annot = False, cmap = 'coolwarm', cbar = False, center=0, ax=ax) #cbar_kws={"orientation": "horizontal"}
ax.set_yticks([])
ax.set_yticklabels([])
ax.set_ylabel('')
ax.set_xticklabels(['Pv','Zm','Sb'],rotation = 360, fontstyle = 'italic')
ax.set_title('-N/Full')

plt.savefig('FN_Heatmap.svg')
plt.savefig('FN_Heatmap.png',dpi=300)
