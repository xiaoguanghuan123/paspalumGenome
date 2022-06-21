import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib import rcParams
import numpy as np
from scipy.stats import ttest_ind

def pcttck(x,pos):
    return "{0}%".format(int((x)))
pter = FuncFormatter(pcttck)

dfcontent = pd.read_csv('DryWeightContent.csv' )
dfcontent = dfcontent.loc[dfcontent['Treatment'] == 'NT']
dfcontent.loc[dfcontent['Condition']=='Full', 'Condition'] = 'Control'
dfcontent.loc[dfcontent['Condition']=='Ndep', 'Condition'] = '-N'
dfcontent.loc[dfcontent['Condition']=='Pdep', 'Condition'] = '-P'
dfN = dfcontent.drop('P', axis = 1)
dfN = dfN


for species, group in dfN.groupby('Species'):
    Full = group.loc[group['Condition'] == 'Control','N']
    NDep= group.loc[group['Condition'] == '-N','N']
    PDep = group.loc[group['Condition'] == '-P','N']
    FullMean = np.mean(Full)
    NDepMean = np.mean(NDep)
    Perc_Red = (FullMean - NDepMean)/FullMean
    tFN = ttest_ind(Full,NDep)
    tFP = ttest_ind(Full,PDep)
    print(species,'N_dep:', tFN[1], Perc_Red)
    print(species,'P_dep:', tFP[1])  
    
    
    
def plotCap(species, vmax, vheight, gap):
    h = vheight + gap
    hh = 2*vheight + gap
    
    xDict = {'Maize':[-0.25, 0, 0.25],'Sorghum':[0.75, 1.0, 1.25],'Paspalum':[1.75, 2.0, 2.25]}
    sigDict = {'Maize':['*', 'ns'],'Sorghum':['**','ns'],'Paspalum':['**','ns']}
    
    x1, x2, x3 = xDict[species][0],xDict[species][1],xDict[species][2]
    t1, t2 = sigDict[species][0], sigDict[species][1]

    ax.plot([x1,x1,x2,x2],[vmax, vmax + vheight, vmax + vheight, vmax], color = 'k', linewidth = 1)
    ax.plot([x1, x1, x3, x3],[vmax + h, vmax + hh, vmax + hh, vmax + h], color = 'k', linewidth = 1)
    
    ax.text((x1 + x2)* .5, vmax + vheight + 0.1, t1, ha ='center', va='center')
    ax.text((x1 + x3)* .5, vmax + hh + 0.1, t2, ha = 'center', va = 'center')
    
    

rcParams['font.family'] = 'sans-serif'
rcParams['font.size'] = 10

myorder = ['Maize','Sorghum','Paspalum']
fig, ax = plt.subplots(figsize=(2.8,3), tight_layout = True)   
plt.setp(ax.spines.values(), linewidth=1.5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', which='major', labelsize=10,direction='out', length=4, width=1.5)

sns.boxplot(x = "Species", y='N', hue="Condition", data=dfN, ax=ax, dodge=True, fliersize=0,linewidth=1, order = myorder)
sns.stripplot(x = "Species", y='N', hue="Condition", data=dfN, dodge=True, linewidth=1, s=2,order = myorder, ax=ax)

handles, labels = ax.get_legend_handles_labels()
labels=['Full','-N', '-P']
l = ax.legend(handles[:3], labels, loc=0, handlelength=0.8, borderpad=0.5)

plt.axvline(0.5,color='black',linestyle= '--', linewidth=1.5)
plt.axvline(1.5,color='black',linestyle= '--', linewidth=1.5)
ax.set_xlim(-0.5,2.5)
ax.set_ylim(-0.5,5)
ax.set_xticklabels(['Zm','Sb','Pv'],size=10,color='k', fontstyle='italic')
ax.set_ylabel('N concentration',size=10, color='black')
ax.set_xlabel('')

#new line of code
ax.yaxis.set_major_formatter(pter)

for asp, grp in dfN.groupby('Species'):
    #print(grp.head())
    vmax = grp['N'].max() + 0.2
    vheight = 0.2
    gap = 0.3
    plotCap(asp, vmax, vheight, gap)

plt.savefig('N_Content.svg')
plt.savefig('N_Content.png', dpi=300)
