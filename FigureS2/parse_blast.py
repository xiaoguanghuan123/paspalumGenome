import os
import pandas as pd

#makeblastdb -in V1_prot.fa -dbtype prot
#blastp -subject Sorgum_bicolor_prot_new.fa -query Pvaginatumv3.1.primaryTrs.pep.fa -evalue 10 -num_threads 8 -outfmt 6 -out PV3SB3.blast

##### BLAST FMT6 HEADER
#query_id        subject_id      pct_identity    aln_length      n_of_mismatches gap_openings    q_start q_end   s_start   s_end   e_value bit_score

header  = open('header.txt').readline().strip().split()
dfbla = pd.read_csv('PV3SB3.blast',sep = '\t',header=None)
dfbla.columns = header
dfbla = dfbla.loc[dfbla['pct_identity']>=60]
print(dfbla.head())


from Bio import SeqIO
def length(fasta):
    lenDict = {}
    fa = SeqIO.parse(fasta, 'fasta')
    for record in fa:
        pid = str(record.id)
        plen = len(str(record.seq))
        lenDict[pid] = plen
    return lenDict

def blastRes(index, dfbla): #'subject_id' 'query_id'
    dfres = pd.DataFrame()
    dfblast = dfbla.set_index(index)
    gList = list(set(dfblast.index.unique()))
    for i, agene in enumerate(gList):
        dfgrp = dfblast.loc[[agene]].reset_index()
        dfgrp = dfgrp.loc[dfgrp.e_value.idxmin()]
        dfres = dfres.append(dfgrp)
    dfres.set_index(index,inplace=True)
    return dfres


pvFa = 'Pvaginatumv3.1.primaryTrs.pep.fa'
lenDict = length(pvFa)
covList = []
dfresPv = blastRes('query_id', dfbla)
for agene in lenDict:
    if not agene in dfresPV.index: continue
    dfgene = dfresPV.loc[[agene]]
    if not len(dfgene.index) == 1: print(agene)
    mLen = float(dfgene['q_end']) - float(dfgene['q_start'])
    coverage = (mLen/lenDict[agene])*100
    covList.append(coverage)
print(len(covList))
print(covList[:20])

sbFa = 'Sorgum_bicolor_prot_new.fa'
lenDict = length(sbFa)
sbcovList = []
dfresSb = blastRes('subject_id', dfbla)

for agene in lenDict:
    if not agene in dfresSb.index: continue
    dfgene = dfresSb.loc[[agene]]
    if not len(dfgene.index) == 1: print(agene)
    mLen = float(dfgene['s_end']) - float(dfgene['s_start'])
    coverage = (mLen/lenDict[agene])*100
    sbcovList.append(coverage)
print(len(sbcovList))
print(sbcovList[:20])   


import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.family'] = 'sans-serif'
rcParams['font.size'] = 10

fig, axes = plt.subplots(2, figsize = (5,5), sharex=True, tight_layout = True)

plt.setp(axes[0].spines.values(), linewidth=1.5)
plt.setp(axes[1].spines.values(), linewidth=1.5)
axes[0].spines['top'].set_visible(False)
axes[0].spines['right'].set_visible(False)
axes[1].spines['top'].set_visible(False)
axes[1].spines['right'].set_visible(False)

axes[0].hist(covList, bins = 30, edgecolor = 'k')
axes[1].hist(sbcovList, bins = 30, edgecolor = 'k')
axes[0].set_ylabel('Number of genes')
axes[1].set_ylabel('Number of genes')
axes[1].set_xlabel('Hit length/gene length')
#ax.set_xticklabels(xticklabels,rotation = 360)
axes[0].set_title('n = {0} (Identity >=60%)'.format(str(len(covList))), fontsize=10)
axes[1].set_title('n = {0} (Identity >=60%)'.format(str(len(sbcovList))),fontsize=10)
plt.savefig('Coverage.svg') 
