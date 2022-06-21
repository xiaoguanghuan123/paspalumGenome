import subprocess as sp
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os, sys


def change_width(ax, new_value) :
    for patch in ax.patches :
        current_width = patch.get_width()
        diff = current_width - new_value
        current_height = patch.get_height()
        # we change the bar width
        patch.set_width(new_value)
        
        # we recenter the bar
        patch.set_x(patch.get_x() + diff * .5)


def add_pvalue(ax, df):
    for i, x in enumerate(plotOrder):
        print(x)
        xpos = i
        ypos = df.loc[df['group']==x,'log2FoldChange'].values * 0.5
        pList = str(df.loc[df['group']==x,'padj'].values[0]).split('e')
        if len(pList) == 2: 
            pVal = 'p = ' + pList[0][:3] + 'e' + pList[-1]
            ax.text(xpos, ypos, pVal, rotation=90, ha ='center', va='bottom')
        else: 
            pVal = 'p = ' + pList[0][:5]
            ax.text(xpos, ypos, pVal, rotation=90, ha ='center', va='bottom')
###############
#trppDict = ['trpp8']

dfNF = pd.read_csv('MaizeValA_NF.csv', index_col = 0)
dfNVF = pd.read_csv('MaizeValA_NVF.csv', index_col = 0)
dfNV = pd.read_csv('MaizeValA_NV.csv', index_col = 0)
dfFV = pd.read_csv('MaizeValA_FV.csv', index_col=0)
dfNF['group'] ='N_VS_F'
dfNVF['group'] = 'NV_VS_F'
dfNV['group'] = 'NV_VS_N'
dfFV['group'] = 'FV_VS_F'




#geneGroup = sys.argv[1] #'trpp11;' #'trehalase' #'trpp6'
#sp.call("grep '{0}' B73v4.gene_function.txt > function_temp.txt".format(geneGroup), shell=True)
dffunc = pd.read_csv("SLC_GENE.txt", sep='\t', header = None, usecols=[0,1])
dffunc.columns = ['geneid','Annotation']
dffunc.set_index('geneid',inplace=True)
dffunc['Annotation'] = dffunc['Annotation'].str.split(';',n=1,expand=True)[0]
print(dffunc)
dfgene = pd.DataFrame()
for agene in dffunc.index:
    if agene in dfNV.index:
        dfgene = dfgene.append(dfNV.loc[agene])
        print('NV_vs_N:' , dffunc.loc[agene]['Annotation'], 'log2FC: ', dfNV.loc[agene]['log2FoldChange'])
for agene in dffunc.index:
    if agene in dfFV.index:
        dfgene = dfgene.append(dfFV.loc[agene])
        print('FV_vs_F:' , dffunc.loc[agene]['Annotation'], 'log2FC: ', dfFV.loc[agene]['log2FoldChange'])
#dfgene.reset_index(inplace=True)
print(dfgene)
plotOrder = ['FV_VS_F','NV_VS_N']
for agene in dfgene.index.unique():
    geneName = dffunc.loc[agene]['Annotation'].split(';')[0].upper().replace('/','_')
    print(agene, geneName)
    df = dfgene.loc[agene]
    print(df)
    fig, ax = plt.subplots(figsize=(1.7,3.4), tight_layout = True) 
    plt.setp(ax.spines.values(), linewidth=1.5)
    plt.setp(ax.spines.values(), linewidth=1.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='both', which='major', labelsize=10, direction='out', length=3, width=1.5)	
    sns.barplot(x = "group", y = 'log2FoldChange',  data=df, order = plotOrder, edgecolor = 'k', linewidth=1, palette = 'Accent', ax=ax)
    change_width(ax, .6)
    add_pvalue(ax, df)        
    ax.set_xticklabels(['Full','-N'])
    ax.set_xlabel('')
    ax.set_ylabel('Log2 fold change (+ValA/-ValA)',size=10)    #$\it{text you want to show in italics}$
    ax.set_title('$\it{0}$'.format(geneName))
    output =  '{0}_fc.png'.format(geneName)
    print(output)
    #plt.show()
    plt.savefig(output)
    #plt.savefig("{0}_fc.png".format(geneName),dpi = 500)


#sns.catplot(kind = 'bar', x = "group", y = 'log2FoldChange',  data=dfgene, order = plotOrder, col = 'index', edgecolor = 'k', dodge=False, linewidth=1, palette = 'Accent', aspect=.8, col_wrap=4)




   
