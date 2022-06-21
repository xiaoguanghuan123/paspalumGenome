import subprocess as sp
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

dfW22 = pd.read_csv('W22_2ndLeaf_NF.csv', index_col = 0)
#dfATG12_1 = pd.read_csv('ATG12-1_2nd_NF.csv', index_col = 0)
dfATG12 = pd.read_csv('ATG12-2_2nd_NF.csv', index_col = 0)
dfW22['group'] ='W22'
dfATG12['group'] = 'ATG12'
print(dfATG12.head())


geneGroup = 'trpp'

sp.call("grep '{0}' B73v4.gene_function.txt > function_temp.txt".format(geneGroup), shell=True)
dffunc = pd.read_csv("function_temp.txt".format(geneGroup), sep='\t', header = None, usecols=[0,1])
dffunc.columns = ['geneid','annotation']
dffunc.set_index('geneid',inplace=True)
dffunc['annotation'] = dffunc['annotation'].str.split(';',n=1,expand=True)[0]
#print(dffunc)
dfgene = pd.DataFrame()
for agene in dffunc.index:
    if agene in dfW22.index:
        dfgene = dfgene.append(dfW22.loc[agene])
        print('W22:' , dffunc.loc[agene]['annotation'], 'log2FC: ', dfW22.loc[agene]['log2FoldChange']) 
for agene in dffunc.index:
    if agene in dfATG12.index:
        dfgene = dfgene.append(dfATG12.loc[agene])
        print('ATG12:' , dffunc.loc[agene]['annotation'], 'log2FC: ', dfATG12.loc[agene]['log2FoldChange'])
dfgene.reset_index(inplace=True)
annList = []
for x in dfgene.index:
    gene = dfgene.loc[x]['index']
    ann = dffunc.loc[gene]['annotation'].strip().upper()
    annList.append(ann)
dfgene['geneName'] = annList


dfheat = dfgene.loc[dfgene['group'] == 'W22', ['log2FoldChange','geneName']]
dfheat['ATG12'] = dfgene.loc[dfgene['group'] == 'ATG12', ['log2FoldChange']].values
dfheat.set_index('geneName',inplace=True)
geneName = list(dfgene.geneName.unique())
geneName.sort(key = lambda x: int(x.strip()[4:]))
#yticklabel = ['TPS' + x[4:] for x in geneName]
yticklabel = ['TPP' + x[4:] for x in geneName]
dfheat = dfheat.loc[geneName]
dfheat.columns = ['W22','ATG12']

fig, ax = plt.subplots(figsize=(2,3.5), tight_layout = True)
plt.setp(ax.spines.values(), linewidth=1.5)
plt.setp(ax.spines.values(), linewidth=1.5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', which='major', labelsize=10, direction='out', length=3, width=1.5)



sns.heatmap(dfheat,annot = False, cmap = 'coolwarm', cbar = False, center=0, ax=ax) #cbar_kws={"orientation": "horizontal"}

##color bar
im = ax.imshow([[-3,3],[-3,3]],cmap=plt.cm.coolwarm)
# create an axes on the right side of ax. The width of cax will be 5%
# of ax and the padding between cax and ax will be fixed at 0.05 inch.
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="10%", pad=0.05)
plt.colorbar(im, cax=cax)
cax.set_title('log2FC',size=8)

ax.set_xticklabels(['WT','ATG12-2'],rotation = 30, ha = 'right')
ax.set_yticklabels(yticklabel, fontstyle='italic',rotation = 360, ha = 'right',fontsize = 10)
ax.set_ylabel('Expression level change (-N/+N)')
plt.savefig('{0}_atg12-2_heatmap.svg'.format(geneGroup))
plt.savefig('{0}_atg12-2_heatmap.png'.format(geneGroup),dpi=500)
