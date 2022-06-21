import pandas as pd
import os 
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles

dfsynList = pd.read_csv('SbZmPvSynListUsedForDEG.csv')
dfsynList.head()

dfSbSyn = dfsynList[['sorghum3','CommonID']]
dfSbSyn.set_index('sorghum3',inplace=True)
dfSbSyn=dfSbSyn[~dfSbSyn.index.duplicated(keep='first')]
SbSynGeneList = set(dfSbSyn.index)
dfSbSyn.info()

dfZmSyn = dfsynList[['CommonID','maize4']]
dfZmSyn.set_index('maize4',inplace=True)
dfZmSyn=dfZmSyn[~dfZmSyn.index.duplicated(keep='first')]
ZmSynGeneList = set(dfZmSyn.index)
dfZmSyn.info()

dfPvSyn =dfsynList[['CommonID','paspalum3']]
dfPvSyn.set_index('paspalum3',inplace=True)
dfPvSyn=dfPvSyn[~dfPvSyn.index.duplicated(keep='first')]
PvSynGeneList = set(dfPvSyn.index)
dfPvSyn.info()

dfZmNF = pd.read_csv('Maize_NF_DEG.csv',sep=',')
dfSbNF = pd.read_csv('Sorghum_NF_DEG.csv',sep=',')
dfPvNF = pd.read_csv('Paspalum_NF_DEG.csv',sep=',')
dfPvNF = dfPvNF.loc[abs(dfPvNF['log2FoldChange'])>=1]
dfZmNF = dfZmNF.loc[abs(dfZmNF['log2FoldChange'])>=1]
dfSbNF = dfSbNF.loc[abs(dfSbNF['log2FoldChange'])>=1]
dfZmNF.rename(columns={dfZmNF.columns[0]:'geneID'},inplace=True)
dfSbNF.rename(columns={dfSbNF.columns[0]:'geneID'},inplace=True)
dfPvNF.rename(columns={dfPvNF.columns[0]:'geneID'},inplace=True)
dfPvNF['geneID'] = dfPvNF.geneID.str.split('.',n=1,expand=True)[0]
print('zm:', len(dfZmNF), 'sb:', len(dfSbNF), 'pv:', len(dfPvNF))

dfZmPF = pd.read_csv('Maize_PF_DEG.csv',sep=',')
dfSbPF = pd.read_csv('Sorghum_PF_DEG.csv',sep=',')
dfPvPF = pd.read_csv('Paspalum_PF_DEG.csv',sep=',')
dfPvPF = dfPvPF.loc[abs(dfPvPF['log2FoldChange'])>=1]

dfZmPF.rename(columns={dfZmPF.columns[0]:'geneID'},inplace=True)
dfSbPF.rename(columns={dfSbPF.columns[0]:'geneID'},inplace=True)
dfPvPF.rename(columns={dfPvPF.columns[0]:'geneID'},inplace=True)
dfPvPF['geneID'] = dfPvPF.geneID.str.split('.',n=1,expand=True)[0]

print('zm:', len(dfZmPF), 'sb:', len(dfSbPF), 'pv:', len(dfPvPF))

ZmNFSyn = [gene for gene in dfZmNF.geneID if gene in ZmSynGeneList]
SbNFSyn = [gene for gene in dfSbNF.geneID if gene in SbSynGeneList]
PvNFSyn = [gene for gene in dfPvNF.geneID if gene in PvSynGeneList]


ZmListNF=set(dfZmSyn.loc[ZmNFSyn].iloc[:,0].values)
SbListNF=set(dfSbSyn.loc[SbNFSyn].iloc[:,0].values)
PvListNF=set(dfPvSyn.loc[PvNFSyn].iloc[:,0].values)


print('ZmNFDEG:',len(dfZmNF.geneID),'ZmNFDEG_syn:',len(ZmNFSyn), 'ZmNFDEG_Nonsyn:',len(dfZmNF.geneID)-len(ZmNFSyn))
print('SbNFDEG:',len(dfSbNF.geneID),'SbNFDEG_syn:',len(SbNFSyn), 'SbNFDEG_Nonsyn:',len(dfSbNF.geneID)-len(SbNFSyn))
print('PvNFDEG:',len(dfPvNF.geneID),'PvNFDEG_syn:',len(PvNFSyn),' PvNFDEG_Nonsyn:',len(dfPvNF.geneID)-len(PvNFSyn))


ZmPFSyn = [gene for gene in dfZmPF.geneID if gene in ZmSynGeneList]
SbPFSyn = [gene for gene in dfSbPF.geneID if gene in SbSynGeneList]
PvPFSyn = [gene for gene in dfPvPF.geneID if gene in PvSynGeneList]


ZmListPF=set(dfZmSyn.loc[ZmPFSyn].iloc[:,0].values)
SbListPF=set(dfSbSyn.loc[SbPFSyn].iloc[:,0].values)
PvListPF=set(dfPvSyn.loc[PvPFSyn].iloc[:,0].values)


print('ZmPFDEG:',len(dfZmPF.geneID),'ZmPFDEG_syn:',len(ZmPFSyn), 'ZmPFDEG_Nonsyn:',len(dfZmPF.geneID)-len(ZmPFSyn))
print('SbPFDEG:',len(dfSbPF.geneID),'SbPFDEG_syn:',len(SbPFSyn), 'SbPFDEG_Nonsyn:',len(dfSbPF.geneID)-len(SbPFSyn))
print('PvPFDEG:',len(dfPvPF.geneID),'PvPFDEG_syn:',len(PvPFSyn),' PvPFDEG_Nonsyn:',len(dfPvPF.geneID)-len(PvPFSyn))


from collections import Counter
dfsynList.set_index('CommonID',inplace=True)
NFDEList =list(ZmListNF)+list(SbListNF)+ list(PvListNF)
DECounter=Counter(NFDEList)
NF3 = [x for x in DECounter if DECounter[x]==3]
PVNF1 = [x for x in PvListNF if DECounter[x] == 1]
SBNF1 = [x for x in SbListNF if DECounter[x] == 1]
ZMNF1 = [x for x in ZmListNF if DECounter[x] == 1]
dfsynList.loc[NF3]['maize4'].to_csv('NF3_ZmID.study',index=False, header = None)
dfsynList.loc[PVNF1]['maize4'].to_csv('PVNF1_ZmID.study',index=False, header = None)
dfsynList.loc[SBNF1]['maize4'].to_csv('SBNF1_ZmID.study',index=False, header = None)
dfsynList.loc[ZMNF1]['maize4'].to_csv('ZMNF1_ZmID.study',index=False, header = None)
print(len(NF3),len(PVNF1))


from collections import Counter
PFDEList=list(ZmListPF)+list(SbListPF)+ list(PvListPF)
PFDECounter=Counter(PFDEList)
PF3 = [x for x in PFDECounter if PFDECounter[x] ==3] 
PVPF1 = [x for x in PvListPF if PFDECounter[x] == 1] 
SBPF1 = [x for x in SbListPF if PFDECounter[x] == 1]
ZMPF1 = [x for x in ZmListNF if PFDECounter[x] == 1]
dfsynList.loc[PF3]['maize4'].to_csv('PF3_ZmID.study',index=False, header = None)
dfsynList.loc[PVPF1]['maize4'].to_csv('PVPF1_ZmID.study',index=False, header = None)
dfsynList.loc[SBPF1]['maize4'].to_csv('SBPF1_ZmID.study',index=False, header = None)
dfsynList.loc[ZMPF1]['maize4'].to_csv('ZMPF1_ZmID.study',index=False, header = None)
print(len(PF3),len(PVPF1))



fpop = open('ZmID_Syn.population','w')
for agene in dfsynList.maize4.unique():
    fpop.write(agene + '\n')
fpop.close()
print('syntenic genes for DEG analysis is:', len(dfsynList.maize4.unique()))


from matplotlib_venn import venn3, venn3_unweighted
fig= plt.figure(figsize=(10,8))
labels = venn3([set(ZmListNF), set(SbListNF), set(PvListNF)],['Maize syntenic DE(-N) genes', 'Sorghum syntenic DE(-N) genes','Paspalum syntenic DE(-N) genes'], alpha=0.5)
venn3_circles([set(ZmListNF), set(SbListNF), set(PvListNF)], linestyle= '-', linewidth=2, color='k')
plt.title('Comparison of nitrogen starvation responsive genes')
plt.savefig('Ndep_SynComparison.svg')


fig= plt.figure(figsize=(10,8))
labels = venn3([set(ZmListPF), set(SbListPF), set(PvListPF)],['Maize syntenic DE(-P) genes', 'Sorghum syntenic DE(-P) genes','Paspalum syntenic DE(-P) genes'], alpha=0.5)
venn3_circles([set(ZmListPF), set(SbListPF), set(PvListPF)], linestyle= '-', linewidth=2, color='k')
plt.title('Comparison of phosphorus starvation responsive genes')
plt.savefig('Pdep_SynComparison.svg')
