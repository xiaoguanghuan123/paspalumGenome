import subprocess as sp
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
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


dffunc = pd.read_csv('SWEET.txt',sep = '\t',index_col=0)
dffunc.head()


dfgene = pd.DataFrame()
for agene in dffunc.index:
    if agene in dfNVF.index:
        dfgene = dfgene.append(dfNVF.loc[agene])
        print('NV_vs_F:' , agene, dffunc.loc[agene]['Annotation'], 'log2FC: ', dfNVF.loc[agene]['log2FoldChange']) 
for agene in dffunc.index:
    if agene in dfNF.index:
        dfgene = dfgene.append(dfNF.loc[agene])
        print('N_vs_F:' ,agene, dffunc.loc[agene]['Annotation'], 'log2FC: ', dfNF.loc[agene]['log2FoldChange'])
for agene in dffunc.index:
    if agene in dfNV.index:
        dfgene = dfgene.append(dfNV.loc[agene])
        print('NV_vs_N:' , agene,dffunc.loc[agene]['Annotation'], 'log2FC: ', dfNV.loc[agene]['log2FoldChange'])
for agene in dffunc.index:
    if agene in dfFV.index:
        dfgene = dfgene.append(dfFV.loc[agene])
        print('FV_vs_F:' , agene, dffunc.loc[agene]['Annotation'], 'log2FC: ', dfFV.loc[agene]['log2FoldChange'])
print(dfgene)

for agene in dffunc.index:
    if not agene in dfgene.index:continue
    geneName = dffunc.loc[agene]['Annotation'].split(';')[0].upper()
    df = dfgene.loc[agene]
    fig, ax = plt.subplots(figsize=(1.7,2.8), tight_layout = True) 
    plt.setp(ax.spines.values(), linewidth=1.5)
    plt.setp(ax.spines.values(), linewidth=1.5)
    #ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='both', which='major', labelsize=10, direction='out', length=3, width=1.5)

    #plotOrder = ['N_VS_F','NV_VS_F','FV_VS_F','NV_VS_N']
    plotOrder = ['FV_VS_F','NV_VS_N']
	
    sns.barplot(x = "group", y = 'log2FoldChange',  data=df, order = plotOrder, edgecolor = 'k', linewidth=1, palette = 'Accent', ax=ax)
    ax.set_xlabel('')
    #ax.set_xticks([0,1,2,3])
    #ax.set_xticklabels(['Control', 'ValA'])
    #ax.set_ylim(df.min()-,1)

    #ax.yaxis.set_major_formatter(pter)

    def change_width(ax, new_value) :
        for patch in ax.patches :
            current_width = patch.get_width()
            diff = current_width - new_value
            current_height = patch.get_height()
            # we change the bar width
            patch.set_width(new_value)

            # we recenter the bar
            patch.set_x(patch.get_x() + diff * .5)
    change_width(ax, .6)


    for i, x in enumerate(plotOrder):
        print(x)
        xpos = i
        ypos = df.loc[df['group']==x,'log2FoldChange'].values - 0.2
        pList = str(df.loc[df['group']==x,'padj'].values[0]).split('e')
        if len(pList) == 2: 
            pVal = 'p = ' + pList[0][:3] + 'e' + pList[-1]
            ax.text(xpos, ypos, pVal, rotation=90, ha ='center', va='bottom')
        else: 
            pVal = 'p = ' + pList[0][:5]
            ax.text(xpos, ypos, pVal, rotation=90, ha ='center', va='bottom')

    ax.set_xticklabels(['Full','-N'])    
    ax.set_ylabel('Expression level change \n log2(+ValA/-ValA)',size=10)    #$\it{text you want to show in italics}$
    ax.set_title('$\it{0}$'.format(geneName)) 

    #ax.set_xticklabels(['Full','-N'], rotation=20, ha='right')    
    #ax.set_ylabel('$\it{0}$ log2 fold change'.format(geneName),size=10)    #$\it{text you want to show in italics}$
    plt.savefig('{0}_fc.svg'.format(geneName))
    plt.savefig("{0}_fc.png".format(geneName),dpi = 500)

