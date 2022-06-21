import pandas as pd
import seaborn as sns
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib import rcParams
import numpy as np

def pcttck(x,pos):
    return "{0}g".format(int(x) )
pter = FuncFormatter(pcttck)

dfbiomass = pd.read_csv('biomass_Fresh_aboveGround.csv', header=None )
dfbiomass.columns = ['Species','Biomass','Condition','Rep']
#dfbiomass.loc[dfbiomass['Condition']=='Full', 'Condition'] = 'Control'

group = dfbiomass.groupby('Species')
for g in group:
    #print(g[1])
    Species = list(set(g[1]['Species']))[0]
    Full = g[1].loc[g[1]['Condition']=='Full', 'Biomass']
    N_dep = g[1].loc[g[1]['Condition']=='-N', 'Biomass']
    P_dep = g[1].loc[g[1]['Condition']=='-P', 'Biomass']
    #print(Full)
    tFN = ttest_ind(Full,N_dep)
    tFP = ttest_ind(Full,P_dep)
    print(Species,'N_dep:', tFN[1])
    print(Species,'P_dep:', tFP[1])  

def plotCap(species, vmax, vheight, gap):
    h = vheight + gap
    hh = 2*vheight + gap
    
    xDict = {'Maize':[-0.25, 0, 0.25],'Sorghum':[0.75, 1.0, 1.25],'Paspalum':[1.75, 2.0, 2.25]}
    sigDict = {'Maize':['***', '***'],'Sorghum':['***','***'],'Paspalum':['ns','ns']}
    
    x1, x2, x3 = xDict[species][0],xDict[species][1],xDict[species][2]
    t1, t2 = sigDict[species][0], sigDict[species][1]

    plt.plot([x1,x1,x2,x2],[vmax, vmax + vheight, vmax + vheight, vmax], color = 'k', linewidth = 1)
    plt.plot([x1, x1, x3, x3],[vmax + h, vmax + hh, vmax + hh, vmax + h], color = 'k', linewidth = 1)
    
    ax.text((x1 + x2)* .5, vmax + vheight + 0.1, t1, ha ='center', va='center')
    ax.text((x1 + x3)* .5, vmax + hh + 0.1, t2, ha = 'center', va = 'center')
    
    

rcParams['font.family'] = 'sans-serif'
rcParams['font.size'] = 10



fig, ax = plt.subplots(figsize=(2.8,3), tight_layout = True)

plt.setp(ax.spines.values(), linewidth=1.5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', which='major', labelsize=10,direction='out', length=4, width=1.5)

sns.boxplot(x="Species", hue="Condition", y="Biomass", data=dfbiomass,ax=ax, dodge=True, fliersize=0,linewidth=1)
sns.stripplot(x="Species", y="Biomass", hue="Condition", data=dfbiomass, dodge=True, linewidth=1, s=2, ax=ax)
#ax.get_legend().remove()
#handles, labels = ax.get_legend_handles_labels()
#l = ax.legend(handles[:3], labels[:3])
handles, labels = ax.get_legend_handles_labels()
labels=['Full','-N', '-P']
l = ax.legend(handles[:3], labels, loc=0, handlelength=0.8, borderpad=0.5)

plt.axvline(0.5,color='black',linestyle= '--', linewidth=1.5)
plt.axvline(1.5,color='black',linestyle= '--', linewidth=1.5)

ax.set_xlim(-0.5,2.5)
ax.set_ylim(-0.5,5.5)

ax.yaxis.set_major_formatter(pter)
ax.set_xticklabels(['Zm','Sb','Pv'],size=10,color='k', fontstyle='italic')
ax.set_ylabel('Above ground fresh biomass',size=10, color='black')
ax.set_xlabel('')


for asp, grp in dfbiomass.groupby('Species'):
    vmax = grp.Biomass.max() + 0.2
    vheight = 0.2
    gap = 0.3
    plotCap(asp, vmax, vheight, gap)

plt.savefig('Biomass_legend.svg')
plt.savefig('Biomass_legend.png', dpi=300)
