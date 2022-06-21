import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib import rcParams
import numpy as np
import scipy 
from scipy import stats
from scipy.stats import ttest_ind


def pcttck(x,pos):
    return "{0}%".format(str(x)[:3])
pter = FuncFormatter(pcttck)


dfcontent = pd.read_csv('DryWeightContent.csv' )
dfcontent = dfcontent.loc[dfcontent['Treatment'] == 'NT']
dfcontent.loc[dfcontent['Condition']=='Full', 'Condition'] = 'Control'
dfcontent.loc[dfcontent['Condition']=='Ndep', 'Condition'] = '-N'
dfcontent.loc[dfcontent['Condition']=='Pdep', 'Condition'] = '-P'
dfP = dfcontent.drop('N', axis = 1)


for species, group in dfP.groupby('Species'):
    Full = group.loc[group['Condition'] == 'Control','P']
    NDep= group.loc[group['Condition'] == '-N','P']
    PDep = group.loc[group['Condition'] == '-P','P']
    FullMean = np.mean(Full)
    PDepMean = np.mean(PDep)
    Perc_Red = (FullMean - PDepMean)/FullMean
    tFN = ttest_ind(Full,NDep)
    tFP = ttest_ind(Full,PDep)
    print(species,'N_dep:', tFN[1])
    print(species,'P_dep:', tFP[1], Perc_Red)  
    
    
def plotCap(species, vmax, vheight, gap):
    h = vheight + gap
    hh = 2*vheight + gap
    
    xDict = {'Maize':[-0.25, 0, 0.25],'Sorghum':[0.75, 1.0, 1.25],'Paspalum':[1.75, 2.0, 2.25]}
    sigDict = {'Maize':['***', '***'],'Sorghum':['***','***'],'Paspalum':['*','*']}
    
    x1, x2, x3 = xDict[species][0],xDict[species][1],xDict[species][2]
    t1, t2 = sigDict[species][0], sigDict[species][1]

    ax.plot([x1,x1,x2,x2],[vmax, vmax + vheight, vmax + vheight, vmax], color = 'k', linewidth = 1)
    ax.plot([x1, x1, x3, x3],[vmax + h, vmax + hh, vmax + hh, vmax + h], color = 'k', linewidth = 1)
    
    ax.text((x1 + x2)* .5, vmax + vheight + 0.01, t1, ha ='center', va='center')
    ax.text((x1 + x3)* .5, vmax + hh + 0.01, t2, ha = 'center', va = 'center')
    
    

rcParams['font.family'] = 'sans-serif'
rcParams['font.size'] = 10

myorder = ['Maize','Sorghum','Paspalum']
fig, ax = plt.subplots(figsize=(3,3), tight_layout = True)   
plt.setp(ax.spines.values(), linewidth=1.5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', which='major', labelsize=10,direction='out', length=4, width=1.5)

sns.boxplot(x = "Species", y='P', hue="Condition", data=dfP, ax=ax, dodge=True, fliersize=0,linewidth=1, order = myorder)
sns.stripplot(x = "Species", y='P', hue="Condition", data=dfP, dodge=True, linewidth=1, s=1.5,order = myorder, ax=ax)

handles, labels = ax.get_legend_handles_labels()
labels=['Full','-N', '-P']
l = ax.legend(handles[:3], labels, loc=0, handlelength=0.8, borderpad=0.5)

plt.axvline(0.5,color='black',linestyle= '--', linewidth=1.5)
plt.axvline(1.5,color='black',linestyle= '--', linewidth=1.5)
ax.set_xlim(-0.5,2.5)
ax.set_ylim(-0.02,0.6)
ax.set_xticklabels(['Zm','Sb','Pv'],size=10,color='k', fontstyle='italic')
ax.set_ylabel('P content per dry weight',size=10, color='black')
ax.set_xlabel('')

ax.yaxis.set_major_formatter(pter)
for asp, grp in dfP.groupby('Species'):
    #print(grp.head())
    vmax = grp['P'].max() + 0.02
    vheight = 0.02
    gap = 0.03
    plotCap(asp, vmax, vheight, gap)

plt.savefig('P_Content.svg')
plt.savefig('P_Content.png', dpi=300)

plt.show()
