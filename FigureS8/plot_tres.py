import os,sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from collections import Counter
import subprocess as sp
###############
#trppDict = ['trpp8']
geneGroup = 'trehalase' #'trehalase' #'trpp6'
sp.call("grep '{0}' B73v4.gene_function.txt > function_{0}.txt".format(geneGroup), shell=True)
dffunc = pd.read_csv("function_{0}.txt".format(geneGroup), sep='\t', header = None, usecols=[0,1])
dffunc.columns = ['geneid','annotation']
dffunc.set_index('geneid',inplace=True)
#dffunc['annotation'] = dffunc['annotation'].str.split(';',n=1,expand=True)[0]
dffunc


synFile = '/home/gsun2unl/Documents/PaspalumGenomeNutrientStressPaper/PaspalumGenome/KaKs/Grass_syntenyList.csv'
spList = ['sorghum3', 'maize1_v4', 'maize2_v4', 'paspalum']
dfsynlist = pd.read_csv(synFile, usecols = spList)
dfsynlist = dfsynlist.loc[dfsynlist['paspalum']!='No Gene']


dfsyn = pd.DataFrame()
for x in dfsynlist.index:
    zm1 = dfsynlist.loc[x]['maize1_v4']
    zm2 = dfsynlist.loc[x]['maize2_v4']
    if zm1 in dffunc.index or zm2 in dffunc.index:
        dfsyn = dfsyn.append(dfsynlist.loc[x])

nameList = []
for x in dfsyn.index:
    geneName = []
    zm1 = dfsyn.loc[x]['maize1_v4']   
    zm2 = dfsyn.loc[x]['maize2_v4']
    if zm1 in dffunc.index: 
        name1 = dffunc.loc[zm1]['annotation'].split(';')[0].strip()
        geneName.append(name1)
    if zm2 != 'No Gene'and zm2 in dffunc.index: 
        name2 = dffunc.loc[zm2]['annotation'].split(';')[0].strip()
        geneName.append(name2)
        print(x, '/'.join(geneName))
    nameList.append('/'.join(geneName))
dfsyn['geneName'] = nameList


dfsyn_melt = dfsyn.melt(id_vars='geneName', var_name='Species', value_name='geneID')
for x in dfsyn_melt.index:
    if not dfsyn_melt.loc[x]['Species'].startswith('maize'):continue
    zmID = dfsyn_melt.loc[x]['geneID']
    if not zmID in dffunc.index: dfsyn_melt.drop(x, inplace=True)
dfsyn_melt.set_index('geneName', inplace=True)
dfsyn_melt.drop_duplicates(inplace=True)
#dfsyn_melt.drop('trpp3',inplace=True)
dfsyn_final = pd.DataFrame()
spDict = {'maize1_v4':'maize1','maize2_v4':'maize2','paspalum':'Paspalum','sorghum3':'Sorghum'}
for x in dfsyn_melt.index.unique():
    dfgene = dfsyn_melt.loc[x]
    varList = dfgene.loc[[x]]['Species'].unique()
    spList = []
    for sp in varList:
        #print(sp)
        spNum = dfgene.Species.value_counts()[sp]
        if  spNum == 1: 
            spList.append(spDict[sp])
        else: 
            spl = [spDict[sp] + str(i+1) for i in range(0,spNum)]
            spList += spl 
    dfsyn_melt.loc[x, 'Species'] = spList
dfsyn_melt

dfZmExp = pd.read_csv('Zm_tpm.csv', sep = '\t', index_col = 0)
dfSbExp = pd.read_csv('Sb_tpm.csv', sep = '\t', index_col = 0)
dfPvExp = pd.read_csv('Pv_tpm.csv', sep = '\t', index_col = 0)
dfZmExp.index=dfZmExp.index.str.rsplit('_', n = 2,expand = True).get_level_values(0)
dfSbExp.index=dfSbExp.index.str.rsplit('.', n = 1,expand = True).get_level_values(0)
dfPvExp.index=dfPvExp.index.str.split('.', n = 1,expand = True).get_level_values(0)

dfexp = pd.DataFrame()
pvList = dfsyn_melt.loc[dfsyn_melt['Species'] == 'Paspalum', 'geneID'].values
zmList = dfsyn_melt.loc[dfsyn_melt['Species'].str.startswith('maize'), 'geneID'].values
sbList = dfsyn_melt.loc[dfsyn_melt['Species'].str.startswith('Sorghum'), 'geneID'].values
dfexp = pd.concat([dfZmExp.loc[zmList], dfSbExp.loc[sbList], dfPvExp.loc[pvList]])
dfexp.columns = dfexp.columns = [x[0] for x in list(dfexp)]
dfsyn_melt = dfsyn_melt.join(dfexp, on='geneID')
dfsyn_melt


for agene in dfsyn_melt.index.unique():
    dfgene = dfsyn_melt.loc[agene].drop('geneID', axis = 1)
    dfgene = dfgene[(dfgene.iloc[:,1:].T >=1).any()]
    plotNum = len(dfgene)
    dfgene = dfgene.melt(id_vars=['Species'], value_name='TPM', var_name = 'condition')
    fig, ax = plt.subplots(figsize=(2.5,3), tight_layout = True) 
    plt.setp(ax.spines.values(), linewidth=1.5)
    plt.setp(ax.spines.values(), linewidth=1.5)
    ax.tick_params(axis='both', which='major', labelsize=10, direction='out', length=3, width=1.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    sns.boxplot(x='Species',y='TPM', hue='condition', order = ['Paspalum', 'maize1', 'Sorghum'], hue_order=['F','N','P'] ,data=dfgene, width = 0.8, fliersize=0, ax=ax)
    sns.stripplot(x='Species',y='TPM', hue='condition', order = ['Paspalum', 'maize1', 'Sorghum'], hue_order=['F','N','P'] ,data=dfgene, s=2, jitter=True, dodge=True, linewidth=1, edgecolor='gray',ax=ax)
    for i in range(0,plotNum):
        xpos = i+0.5 
        ax.axvline(xpos,color='black',linestyle= '--', linewidth=1)
    interval = np.linspace(0, dfgene.TPM.max()+5, 5, dtype=int)
    #print(interval)
    
    ax.set_ylim(0, dfgene.TPM.max() + interval[1])
    handles, labels = ax.get_legend_handles_labels()
    labels=['Full','-N','-P']
    l = ax.legend(handles[:3], labels, loc=0, handlelength=0.8, borderpad=0.5)
    ax.set_xlabel('')
    ax.set_xticklabels(['Pv','Zm', 'Sb'], fontstyle = 'italic')
    ax.set_title('trehalase')
    ax.set_ylabel('TPM',size=10, color='black')
    plt.savefig('{0}_tpm.svg'.format(agene.replace('/','_')))
