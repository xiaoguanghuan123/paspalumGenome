import subprocess as sp
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

dfFV = pd.read_csv('MaizeValA_FV.csv', index_col = 0)
dfNV = pd.read_csv('MaizeValA_NV.csv', index_col = 0)
dfFV['group'] ='FV'
dfNV['group'] = 'NV'



geneGroup = 'SLC'
#sp.call("grep '{0}' B73v4.gene_function.txt > function_temp.txt".format(geneGroup), shell=True)
dffunc = pd.read_csv("SLC_GENE.txt", sep='\t', header = None, usecols=[0,1])
dffunc.columns = ['geneid','Annotation']
dffunc.set_index('geneid',inplace=True)
dffunc['Annotation'] = dffunc['Annotation'].str.split(';',n=1,expand=True)[0]
print(dffunc)
dfgene = pd.DataFrame()
for agene in dffunc.index:
    if agene in dfFV.index:
        dfgene = dfgene.append(dfFV.loc[agene])
        #print('W22:' , dffunc.loc[agene]['annotation'], 'log2FC: ', dfW22.loc[agene]['log2FoldChange']) 
for agene in dffunc.index:
    if agene in dfNV.index:
        dfgene = dfgene.append(dfNV.loc[agene])
        #print('ATG12:' , dffunc.loc[agene]['annotation'], 'log2FC: ', dfATG12.loc[agene]['log2FoldChange'])
dfgene.reset_index(inplace=True)
annList = []
for x in dfgene.index:
    gene = dfgene.loc[x]['index']
    ann = dffunc.loc[gene]['Annotation'].strip()
    annList.append(ann)
dfgene['geneName'] = annList


dfheat = dfgene.loc[dfgene['group'] == 'FV', ['log2FoldChange','geneName']]
dfheat['NV'] = dfgene.loc[dfgene['group'] == 'NV', ['log2FoldChange']].values
dfheat.set_index('geneName',inplace=True)
geneName = list(dfgene.geneName.unique())
geneName.sort()
print(geneName)
#geneName.sort(key = lambda x: int(x.strip()[4:]))
#yticklabel = ['TPS' + x[4:] for x in geneName]
#yticklabel = ['ATG' + x[4:] for x in geneName]
yticklabel = geneName
yticklabel
dfheat = dfheat.loc[geneName]
dfheat.columns = ['Full','-N']

fig, ax = plt.subplots(figsize=(2.5,3.5), tight_layout = True)
plt.setp(ax.spines.values(), linewidth=1.5)
plt.setp(ax.spines.values(), linewidth=1.5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', which='major', labelsize=10, direction='out', length=3, width=1.5)



sns.heatmap(dfheat,annot = False, cmap = 'coolwarm', cbar = False, center=0, ax=ax) #cbar_kws={"orientation": "horizontal"}

'''
##color bar
im = ax.imshow([[-1.5,1.5],[-1.5,1.5]],cmap=plt.cm.coolwarm)
# create an axes on the right side of ax. The width of cax will be 5%
# of ax and the padding between cax and ax will be fixed at 0.05 inch.
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="20%", pad=0.1)
plt.colorbar(im, cax=cax)
cax.set_title('Log2FC',size=10)
'''
ax.set_xticklabels(['Full','-N'],rotation = 30, ha = 'right')
ax.set_yticklabels(yticklabel, fontstyle='italic',rotation = 360, ha = 'right',fontsize = 10)
ax.set_ylabel('Log2 fold change (+ValA/-ValA)')
plt.savefig('{0}_B73_heatmap.svg'.format(geneGroup))
plt.savefig('{0}_B73_heatmap.png'.format(geneGroup),dpi=500)
