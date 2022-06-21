import subprocess as sp
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rcParams
import os

fcDir = 'FoldChange_ttest_GS'
fcList = [afile for afile in os.listdir(fcDir) if afile.startswith('Combined')]
fcNList = [os.path.join(fcDir, afile) for afile in fcList if 'Ndep' in afile]
#fcPList = [os.path.join(fcDir, afile) for afile in fcList if 'Pdep' in afile]

dffc_N = pd.DataFrame()
for i, afile in enumerate(fcNList):
    print(afile)
    sp = afile.split('/')[1].split('_')[2]
    dfNfc = pd.read_csv(afile, sep = '\t', index_col = 0)
    dfNfc = dfNfc.iloc[:,[1,2,4]]
    #print(dfNfc.head())
    foldChange = []
    for x in dfNfc.index:
        if dfNfc.loc[x]['p.value']>0.05:foldChange.append(0)
        elif abs(dfNfc.loc[x]['log2(FC)']) < 1 : foldChange.append(0)
        else: foldChange.append(dfNfc.loc[x]['log2(FC)'])
    dfNfc[sp] = foldChange
    if i == 0: 
        dffc_N = pd.DataFrame(dfNfc[sp])
    else:
        dfNfc = dfNfc.loc[dffc_N.index]
        dffc_N[sp] = dfNfc[sp]
        
        
for x in dffc_N.index:
    i=0
    for v in dffc_N.loc[x].values:
        if abs(v) <= 1: i += 1
    if i == 3: dffc_N = dffc_N.drop(x, axis = 0)

nsList = []
for x in dffc_N.index:
    vCount = dict(dffc_N.loc[x].value_counts())
    if 0 in vCount: nsList.append(vCount[0])
    else: nsList.append(0)

dffc_N['nsCount'] = nsList
dffc_N.sort_values(by = 'nsCount', ascending = True, inplace=True)
dffc_N = dffc_N.drop('nsCount',axis =1)
dffc_N

from matplotlib import rcParams
fig, ax = plt.subplots(figsize=(12,3),tight_layout = True)
sns.heatmap(dffc_N.transpose(), annot = False, cmap = 'coolwarm', cbar = False, center=0, ax=ax) #cbar_kws={"orientation": "horizontal"}
rcParams['font.family'] = 'sans-serif'
rcParams['font.size'] = 10
ax.set_xticklabels(dffc_N.index, rotation = 35, ha = 'right')
ax.set_yticklabels(['P. vaginatum','Z. mays', 'S. bicolor'], rotation = 360, fontstyle = 'italic')

# create an axes on the right side of ax. The width of cax will be 5%
# of ax and the padding between cax and ax will be fixed at 0.05 inch.
'''
im = ax.imshow([[-3,3],[-3,3]],cmap=plt.cm.coolwarm)
divider = make_axes_locatable(ax)
cax = divider.append_axes("top", size="20%", pad="20%",)
fig.colorbar(im, orientation="horizontal", cax=cax)
#cax.set_title('Fold change',size=8)
'''
plt.savefig('Ndep_heatmap.svg')
plt.savefig('Ndep_heatmap.png',dpi = 500)
