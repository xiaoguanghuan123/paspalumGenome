import subprocess as sp
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

dfNF = pd.read_csv('MaizeValA_NF.csv', index_col = 0)
dfNVF = pd.read_csv('MaizeValA_NVF.csv', index_col = 0)
dfNV = pd.read_csv('MaizeValA_NV.csv', index_col = 0)
dfFV = pd.read_csv('MaizeValA_FV.csv', index_col=0)
dfNF['group'] ='N_VS_F'
dfNVF['group'] = 'NV_VS_F'
dfNV['group'] = 'NV_VS_N'
dfFV['group'] = 'FV_VS_F'


dffunc = pd.read_csv('gln_gene_models.txt',sep = '\t', usecols=['gene_id', 'locus_name'])
dffunc.columns = ['geneid','annotation']
dffunc.set_index('geneid',inplace=True)
print(dffunc.head())

dfgene = pd.DataFrame()
for agene in dffunc.index:
    if agene in dfFV.index:
        dfgene = dfgene.append(dfFV.loc[agene])
        print('FV:' , dffunc.loc[agene]['annotation'], 'log2FC: ', dfFV.loc[agene]['log2FoldChange']) 
for agene in dffunc.index:
    if agene in dfNV.index:
        dfgene = dfgene.append(dfNV.loc[agene])
        print('NV:' , dffunc.loc[agene]['annotation'], 'log2FC: ', dfNV.loc[agene]['log2FoldChange'])
dfgene.reset_index(inplace=True)
annList = []
for x in dfgene.index:
    gene = dfgene.loc[x]['index']
    ann = dffunc.loc[gene]['annotation'].strip().upper()
    annList.append(ann)
dfgene['geneName'] = annList


dfheat = dfgene.loc[dfgene['group'] == 'FV_VS_F', ['log2FoldChange','geneName']]
dfheat['NV'] = dfgene.loc[dfgene['group'] == 'NV_VS_N', ['log2FoldChange']].values
dfheat.set_index('geneName',inplace=True)
geneName = list(dfgene.geneName.unique())
geneName.sort(key = lambda x: int(x.strip()[3:]))
print(geneName)
yticklabel = ['GLN' + x[3:] for x in geneName]
dfheat = dfheat.loc[geneName]
dfheat.columns = ['Full','-N']

fig, ax = plt.subplots(figsize=(2,3.5), tight_layout = True)
plt.setp(ax.spines.values(), linewidth=1.5)
plt.setp(ax.spines.values(), linewidth=1.5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', which='major', labelsize=10, direction='out', length=3, width=1.5)


sns.heatmap(dfheat,annot = False, cmap = 'coolwarm', cbar = False, center=0, ax=ax) #cbar_kws={"orientation": "horizontal"}

##color bar
im = ax.imshow([[-3,3],[-3,3]],cmap=plt.cm.coolwarm)
# create an axes on the right side of ax. The width of cax will be 5% of ax and the padding between cax and ax will be fixed at 0.05 inch.
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="15%", pad="10%")
plt.colorbar(im, cax=cax)
cax.set_title('Log2FC', size = 10)

ax.set_xticklabels(['Full','-N'],rotation = 30, ha = 'right')
ax.set_yticklabels(yticklabel, fontstyle='italic',rotation = 360, ha = 'right',fontsize = 10)
ax.set_ylabel('Log2 fold change (+ValA/-ValA)')
plt.savefig('Gln_ValA_heatmap.svg')
plt.savefig('Gln_ValA_heatmap.png', dpi=500)



