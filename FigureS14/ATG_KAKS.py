import os 
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import subprocess as sp



dfkaks = pd.read_csv('SevenSpeciesKAKS_V3.csv',sep='\t', index_col = 0)
dfSynMerge = pd.read_csv('Syntenic_GENE_SET.csv', index_col = 0)


geneGroup ='autophagy'
geneGroup2 ='Autophagy'
sp.call("grep '{0}' B73v4.gene_function.txt > function_{1}.txt".format(geneGroup, geneGroup), shell=True)
sp.call("grep '{0}' B73v4.gene_function.txt > function_{1}.txt".format(geneGroup2, geneGroup2), shell=True)
sp.call("cat function_{0}.txt function_{1}.txt > function_autophagy_combined.txt".format(geneGroup, geneGroup2), shell=True)
dfTrefunc = pd.read_csv('function_autophagy_combined.txt'.format(geneGroup), sep='\t', header = None, usecols=[0,1])
dfTrefunc.columns = ['geneid','annotation']
dfTrefunc.set_index('geneid',inplace=True)
dfTrefunc['annotation'] = dfTrefunc['annotation'].str.split(';',n=1,expand=True)[0]

trpsList = list(dfTrefunc.index)
trpsList = [x for x in dfSynMerge.index if dfSynMerge.loc[x]['Maize'] in trpsList] 
geneList = dfSynMerge.loc[trpsList, 'Maize']

dftreSynkaks = dfkaks.loc[trpsList]
dftreSynkaks['B73v4_geneID'] = geneList
dftreSynkaks.to_csv('{0}GeneOrtho_kaks.csv'.format(geneGroup))



import matplotlib as mpl
import matplotlib.font_manager as font_manager
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.size'] = 10

from scipy.stats import wilcoxon
fig, ax = plt.subplots(figsize=(2.7 ,3.5),tight_layout = True)
plt.setp(ax.spines.values(), linewidth=1.5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', which='major', labelsize=10,direction='out', length=3, width=1.5)

l1 = ['Pavag','Seita', 'Sobic','Ortho']

for i,species in enumerate(l1[1:]):
    Pv = dftreSynkaks['Pavag'].values
    Sp = dftreSynkaks[species].values
    w, p = wilcoxon(Pv,Sp)
    print(species, w, p)
    
    
dfmelt = dftreSynkaks[l1].melt()
dfmelt.columns = ['Species', 'kaks']
dfmelt.to_csv('{0}SynKAKS_MELTED.csv'.format(geneGroup))
print(dfmelt.head())
sns.boxplot(x = "Species", y = "kaks",  dodge=False, data=dfmelt,fliersize=0, boxprops=dict(alpha = 0.5), width = 0.3, linewidth=1, ax=ax)
sns.stripplot(x = "Species", y = "kaks",data=dfmelt, s=3, ax=ax)

ax.set_ylim(-0.02,0.65)
ax.set_ylabel('Ka/Ks Ratio',size=10)
ax.set_title('n=23')#[:3]+'e'+str(p).split('e')[1])
ax.set_xticklabels(['Pv', 'Si', 'Sb', 'Ot'],size=10,color='k', fontstyle='italic') #, rotation=20, ha = 'right')
ax.set_xlabel('')


#plot horizontal lines for the bars
plt.plot([0,1], [0.5, 0.5], linewidth=1.5, solid_joinstyle = 'miter', color='k')
plt.plot([0,2], [0.55, 0.55], linewidth=1.5, solid_joinstyle = 'miter', color='k')
plt.plot([0,3], [0.6, 0.6], linewidth=1.5, solid_joinstyle = 'miter', color='k')


#plot vertical ticks for the bars
plt.plot([0,0], [0.49, 0.50], linewidth=1.5, solid_joinstyle = 'miter', color='k')
plt.plot([0,0], [0.54, 0.55], linewidth=1.5, solid_joinstyle = 'miter', color='k')
plt.plot([0,0], [0.59, 0.60], linewidth=1.5, solid_joinstyle = 'miter', color='k')

plt.plot([1,1], [0.49, 0.50], linewidth=1.5, solid_joinstyle = 'miter', color='k')
plt.plot([2,2], [0.54, 0.55], linewidth=1.5, solid_joinstyle = 'miter', color='k')
plt.plot([3,3], [0.59, 0.60], linewidth=1.5, solid_joinstyle = 'miter', color='k')

#plt.text(s='p = 0.55', x=0.4, y = 0.51)
plt.savefig('{0}GeneKAKS.svg'.format(geneGroup))
plt.savefig('{0}GeneKAKS.png'.format(geneGroup), dpi = 500)
