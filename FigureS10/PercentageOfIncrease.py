import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib import rcParams
import seaborn as sns
rcParams['font.family'] = 'sans-serif'
rcParams['font.size'] = 10


def pcttck(x,pos):
    return "{0}%".format(str(x)[:3])
pter = FuncFormatter(pcttck)

nDec, pDec = [], []
spList = ['Maize','Sorghum','Paspalum']
spList = ['Maize','Paspalum']
ind = [0,1]
width = 1

dfcontent = pd.read_csv('DryWeightContent.csv')
#dfcontent = dfcontent.loc[dfcontent['Treatment'] == 'NT']
dfcontent.loc[dfcontent['Condition']=='Full', 'Condition'] = 'Control'
dfcontent.loc[dfcontent['Condition']=='Ndep', 'Condition'] = '-N'
dfcontent.loc[dfcontent['Condition']=='Pdep', 'Condition'] = '-P'
dfN = dfcontent.drop('P', axis = 1)
dfP = dfcontent.drop('N', axis = 1)

print(dfN.head())

dfspP = dfP.loc[(dfP['Species'] != 'Sorghum') & (dfP['Condition'] == 'Control')]
dfspN = dfN.loc[(dfN['Species'] != 'Sorghum') & (dfN['Condition'] == 'Control')]

import scipy
from scipy import stats

for species, grp in dfspP.groupby('Species'):
	nt = grp.loc[grp['Treatment']=='NT', 'P'].values
	valA = grp.loc[grp['Treatment']=='ValA', 'P'].values
	t, p = scipy.stats.ttest_ind(nt, valA)
	print(species, 'P content', 'Pval:', p)
for species, grp in dfspN.groupby('Species'):
	nt = grp.loc[grp['Treatment']=='NT', 'N'].values
	valA = grp.loc[grp['Treatment']=='ValA', 'N'].values
	t, p = scipy.stats.ttest_ind(nt, valA)
	print(species, 'N content', 'Pval:', p)


fig, axes = plt.subplots(1,2, figsize=(4,2.9), tight_layout = True)   
plt.setp(axes[0].spines.values(), linewidth=1.5)
plt.setp(axes[1].spines.values(), linewidth=1.5)
axes[0].spines['top'].set_visible(False)
axes[0].spines['right'].set_visible(False)
axes[0].tick_params(axis='both', which='major', labelsize=10, direction='out', length=4, width=2)
axes[1].spines['top'].set_visible(False)
axes[1].spines['right'].set_visible(False)
axes[1].tick_params(axis='both', which='major', labelsize=10, direction='out', length=4, width=2)


sns.boxplot(x = "Species", y='P', hue="Treatment", data=dfspP, dodge=True, boxprops=dict(alpha = 0.5), fliersize=0, linewidth=1, palette = 'Set2', order = spList, ax=axes[0])
sns.stripplot(x = "Species", y='P', hue="Treatment", data=dfspP, dodge = True,  jitter = 0.02,  s=3,  order = spList, palette = 'Set2', ax=axes[0])

sns.boxplot(x = "Species", y='N', hue="Treatment", data=dfspN, dodge=True, boxprops=dict(alpha = 0.5), fliersize=0, linewidth=1, palette = 'Set2', order = spList, ax=axes[1])
sns.stripplot(x = "Species", y='N', hue="Treatment", data=dfspN, dodge = True,  jitter = 0.02,  s=3,  order = spList, palette = 'Set2', ax=axes[1])


axes[0].get_legend().remove()
axes[1].get_legend().remove()

axes[0].set_xlabel('')
axes[1].set_xlabel('')

axes[0].set_ylabel('P content per dry weight')
axes[1].set_ylabel('N content per dry weight')

axes[1].set_ylim(1.5,2.6)
axes[0].yaxis.set_major_formatter(pter)
axes[1].yaxis.set_major_formatter(pter)

axes[0].set_xticklabels(['Z. mays', 'P. vaginatum'], rotation = 20, ha = 'right', fontstyle = 'italic')
axes[1].set_xticklabels(['Z. mays', 'P. vaginatum'], rotation = 20, ha = 'right', fontstyle = 'italic')


def plotCap(species, vmax, vheight, gap, t1, ax):
    h = vheight + gap
    xDict = {'Maize':[-0.2, 0.2],'Paspalum':[0.8, 1.2]}
    x1, x2 = xDict[species][0],xDict[species][1]
    plt.plot([x1,x1,x2,x2],[vmax, vmax + vheight, vmax + vheight, vmax], color = 'k', linewidth = 1)
    ax.text((x1 + x2)* .5, vmax + vheight + 0.01, t1, ha ='center', va='center')
 
sigDict = {'Maize':{'P':'***','N':'*'},'Paspalum':{'P':'ns','N':'ns'}}   
xDict = {'Maize':[-0.2, 0.2],'Paspalum':[0.8, 1.2]}
for species, grp in dfspP.groupby('Species'):
    print(grp)
    vmax = grp['P'].max() + 0.005
    vheight = 0.002
    gap = 0.001
    t1 = sigDict[species]['P']
    ax = axes[0]
    h = vheight + gap
    
    x1, x2 = xDict[species][0],xDict[species][1]
    ax.plot([x1,x1,x2,x2],[vmax, vmax + vheight, vmax + vheight, vmax], color = 'k', linewidth = 1)
    ax.text((x1 + x2)* .5, vmax + vheight + 0.01, t1, ha ='center', va='center')

for species, grp in dfspN.groupby('Species'):
    print(grp)
    vmax = grp['N'].max() + 0.025
    vheight = 0.02
    gap = 0.015
    t1 = sigDict[species]['N']
    ax = axes[1]
    h = vheight + gap
    
    x1, x2 = xDict[species][0],xDict[species][1]
    ax.plot([x1,x1,x2,x2],[vmax, vmax + vheight, vmax + vheight, vmax], color = 'k', linewidth = 1)
    ax.text((x1 + x2)* .5, vmax + vheight + 0.05, t1, ha ='center', va='center')



'''
for species in spList:
    dfspP = dfP.loc[(dfP['Species'] == species) & (dfP['Condition'] == 'Control')]
    dfspN = dfN.loc[(dfN['Species'] == species) & (dfN['Condition'] == 'Control')]
    print(dfspN)
    NT_Pi = dfspP.loc[dfspP['Treatment'] == 'NT','P'].mean()
    NT_Ni = dfspN.loc[dfspN['Treatment'] == 'NT','N'].mean()


    ValA_Pi = dfspP.loc[dfspP['Treatment'] == 'ValA','P'].mean()
    ValA_Ni = dfspN.loc[dfspN['Treatment'] == 'ValA','N'].mean()
    

    Pi_Perc_Red = 100*(ValA_Pi - NT_Pi)/ValA_Pi
    Ni_Perc_Red = 100*(ValA_Ni - NT_Ni)/ValA_Ni
    
    #print(species, 'Pi:', Pi_Perc_Red)
    #print(species, 'Ni:', Ni_Perc_Red)
    nDec.append(Ni_Perc_Red)
    pDec.append(Pi_Perc_Red)

    

axes[0].bar(ind, nDec, width = 0.5, edgecolor = 'k', color = 'cyan')
#axes[0].invert_yaxis()
axes[1].bar(ind, pDec, width = 0.5, edgecolor = 'k', color = 'cyan')
#axes[1].invert_yaxis()
axes[0].yaxis.set_major_formatter(pter)
axes[1].yaxis.set_major_formatter(pter)

axes[0].set_xticks(ind)
axes[0].set_xticklabels(['Z. mays','S. bicolor','P. vaginatum'],size=10,color='k', fontstyle='italic', rotation=20, ha = 'right')
axes[1].set_xticks(ind)
axes[1].set_xticklabels(['Z. mays','S. bicolor','P. vaginatum'],size=10,color='k', fontstyle='italic', rotation=20, ha = 'right')

axes[0].set_ylabel('Increase in N')
axes[1].set_ylabel('Increase in P')
'''
plt.savefig('PercentageOfIncrease.svg')
plt.savefig("PercentageOfIncrease.png",dpi = 500)
