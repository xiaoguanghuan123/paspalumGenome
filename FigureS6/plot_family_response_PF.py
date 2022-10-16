import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

OrthoFinderDir = 'Orthogroups.txt'
dfFamilyMCL = pd.read_csv('Orthogroups.txt', sep = ':', header = None, index_col = 0)
dfFamilyMCL.index = list(map(int,(dfFamilyMCL.index.str[4:])))
dfFamilyMCL.index = dfFamilyMCL.index + 1
fDict_PV = {}
for afamily in dfFamilyMCL.index:
    #familyID = 'f' + str(afamily)
    gList = dfFamilyMCL.loc[afamily][1].split()
    PvGeneList = [x[5:] for x in gList if x.startswith('Pvag')]
    if not len(PvGeneList)==0: fDict_PV[afamily] = PvGeneList
dfPvFamily = pd.DataFrame.from_dict(fDict_PV,orient = 'index')
dfPvFamily

expandFamilyStudy = open('FamilyExpandedInPV.study').readlines()
expandFamilyList=[int(x.strip()) for x in expandFamilyStudy]

PvFunction = open('gene.functions.txt','r')
aDict={}
for aline in PvFunction:
    x = aline.strip().split('\t')
    geneID = x[0].strip('.g')
    if x[2] == 'PFAM': 
        function = x[-1]
        if not geneID in aDict: aDict[geneID]=[function]
        else: aDict[geneID].append(function)
        
        
dfPvExpandGene = dfPvFamily.loc[expandFamilyList]
aList=[]
gList=[]
for x in dfPvExpandGene.index:
    geneList = list(dfPvExpandGene.loc[x].dropna().values)
    faList =[','.join(aDict[geneID]) for geneID in geneList if geneID in aDict]
    faList = list(set(faList))
    aList.append(';'.join(faList))
    g1List = ';'.join(geneList)
    gList.append(g1List)
dfPvExpandGene['function']=aList
dfPvExpandGene['geneid'] = gList
dfPvExpandGene = dfPvExpandGene[['function','geneid']]
#dfPvExpandGene.reset_index()['index'].to_csv(goatoolsDir + 'GeneFamilyDEG.pop',header=None, index=False)
dfPvExpandGene

fList=[]
for a in dfPvExpandGene.function:
    aList = [x.lower() for x in a.split(';')]
    aList = set(aList)
    function = ';'.join(aList)
    fList.append(function)
dfPvExpandGene['function']=fList


PvDEGNF = 'Paspalum_NF_DEG.csv'
PvDEGPF = 'Paspalum_PF_DEG.csv'
dfPvDEGNF = pd.read_csv(PvDEGNF,sep=',' , index_col = 0)
dfPvDEGPF = pd.read_csv(PvDEGPF,sep=',' , index_col = 0)
dfPvDEGNF.head()

#extract the gene families with genes that are responding to nutrient starvation
NFDEG=set(dfPvDEGNF.index)
PFDEG=set(dfPvDEGPF.index)
PvExpandGeneDEG={}
PvExpandGeneDEG['Ndep']={}
PvExpandGeneDEG['Pdep']={}
for afamily in dfPvExpandGene.index:
    PvExpandGeneDEG['Ndep'][afamily]=[]
    PvExpandGeneDEG['Pdep'][afamily]=[]
    gList = dfPvExpandGene.loc[afamily].geneid.split(';')
    NFList = [agene + '.g' for agene in gList if agene + '.g' in NFDEG]
    PFList = [agene + '.g' for agene in gList if agene + '.g' in PFDEG]
    if len(NFList) !=0: PvExpandGeneDEG['Ndep'][afamily].append(';'.join(NFList))
    if len(PFList) !=0: PvExpandGeneDEG['Pdep'][afamily].append(';'.join(PFList))
dfPvExpandNFDEG = pd.DataFrame.from_dict(PvExpandGeneDEG['Ndep'],orient='index',columns=['geneid']).dropna()
dfPvExpandPFDEG = pd.DataFrame.from_dict(PvExpandGeneDEG['Pdep'],orient='index',columns=['geneid']).dropna()
dfPvExpandNFDEG['function']=dfPvExpandGene.loc[dfPvExpandNFDEG.index]['function']
dfPvExpandPFDEG['function']=dfPvExpandGene.loc[dfPvExpandPFDEG.index]['function']



dfPvtpm = pd.read_csv('Pv_tpm.csv',sep='\t',index_col = 0)
dfPvtpm.index= [agene + '.g' for agene in dfPvtpm.index] 
dfPvtpm

def meltdataframeForboxplot(fid,df):
    geneList = df.loc[fid].geneid.split(';')
    familyid = df.loc[fid].name.astype(str)
    df = dfPvtpm.loc[geneList].transpose()
    for a in list(df):
        if df[a].mean()<1.5:df.drop(labels=a,axis=1,inplace=True)
    dfmelt = df.reset_index().melt(id_vars=['index'])
    v=dfmelt['index']
    condition = [x[0] for x in v]
    dfmelt['condition']=condition
    return dfmelt
dfout = pd.DataFrame()
fig,axes = plt.subplots(1, 5, figsize=(14,3.5), tight_layout =True)
a=0
for afamily in dfPvExpandPFDEG.index: 
    function = dfPvExpandGene.loc[afamily]['function']
    #if 'domain of unknown function' in function: continue
    #if function == '':continue
    #if afamily == 339: continue
    
    dfmelt = meltdataframeForboxplot(afamily,dfPvExpandPFDEG)
    
    if not len(dfmelt)==0:
        
        dfgroupMean = dfmelt.groupby('condition').mean()  
        Nfpkm = dfgroupMean.loc['N'].value
        Ffpkm = dfgroupMean.loc['F'].value
        if Nfpkm < 1 :continue
        elif Nfpkm/Ffpkm < 2: continue
        #row = divmod(int(a),5)[0]
        col = divmod(int(a),5)[1]
        a+=1
        #print(row,col)
        #dfmelt = dfmelt.loc[dfmelt.condition!='P']
        dfmelt.rename(columns={'variable':'Gene_id'},inplace=True)
        ax = axes[col]
        sns.pointplot(x="condition", y="value", hue="Gene_id", data=dfmelt, linewidth=3,dodge=True, ax=ax)
        plt.setp(ax.spines.values(), linewidth=1.5)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(axis='both', which='major', labelsize=10, direction='out', length=3, width=1.5)
        #ax.set_xlabel('Condition',size=10)
        ax.set_ylabel('TPM',size=10)
        title = function.split(';')[0]
        ax.set_xticklabels(['Full','-N','-P'])
        
        #dfmelt.boxplot(by='condition')
        if afamily==7:
            title = 'Leucine rich repeat kinase'
            ax.set_title(title , fontsize=10)
        elif afamily==71:
            title = 'wall-associated receptor kinase'
            ax.set_title(title , fontsize=10)
        elif afamily==144:
            title = 'NADH oxidase family'
            ax.set_title(title , fontsize=10)
        elif afamily==204:
            title = 'Leucin-rich kinase domain;'
            ax.set_title(title , fontsize=10)
        elif afamily==272:
            title = 'Terpene synthase'
            ax.set_title(title , fontsize=10)
        else:
            title = function.split(';')[0]
            ax.set_title(title , fontsize=10)
            print(afamily, title)
            
        dfout = dfout.append(dfmelt)
        #print(afamily,dfPvExpandGene.loc[afamily]['function'])

dfout.to_csv('PF_FamilyResponse.csv')
plt.savefig('PvSpecificExpandedFamilyPFDEG.svg')
plt.savefig('PvSpecificExpandedFamilyPFDEG.png',dpi=300)
