import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib import rcParams
import numpy as np

def pcttck(x,pos):
    return "{0}cm".format(int(x) )
pter = FuncFormatter(pcttck)

dfbiomass = pd.read_csv('RootLength.csv' )
dfbiomass.loc[dfbiomass['Condition']=='Full', 'Condition'] = 'Control'
print(dfbiomass.head())

def plotCap(species, vmax, vheight, gap):
    h = vheight + gap
    hh = 2*vheight + gap
    
    xDict = {'Maize':[-0.25, 0, 0.25],'Sorghum':[0.75, 1.0, 1.25],'Paspalum':[1.75, 2.0, 2.25]}
    sigDict = {'Maize':['****', 'ns'],'Sorghum':['***','*'],'Paspalum':['*','ns']}
    
    x1, x2, x3 = xDict[species][0],xDict[species][1],xDict[species][2]
    t1, t2 = sigDict[species][0], sigDict[species][1]

    plt.plot([x1,x1,x2,x2],[vmax, vmax + vheight, vmax + vheight, vmax], color = 'k', linewidth = 1)
    plt.plot([x1, x1, x3, x3],[vmax + h, vmax + hh, vmax + hh, vmax + h], color = 'k', linewidth = 1)
    
    ax.text((x1 + x2)* .5, vmax + vheight + 0.1, t1, ha ='center', va='bottom')
    ax.text((x1 + x3)* .5, vmax + hh + 0.1, t2, ha = 'center', va = 'bottom')
    
    

rcParams['font.family'] = 'sans-serif'
rcParams['font.size'] = 10

fig, ax = plt.subplots(figsize=(3.1,3), tight_layout = True)
plt.setp(ax.spines.values(), linewidth=1.5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', which='major', labelsize=10,direction='out', length=4, width=1.5)

sns.boxplot(x="Species", hue="Condition", y="Length", data=dfbiomass, ax=ax, dodge=True, fliersize=0, linewidth=1)
sns.stripplot(x="Species", y="Length", hue="Condition", data=dfbiomass, dodge=True, linewidth=1, s=2, ax=ax)
ax.get_legend().remove()
handles, labels = ax.get_legend_handles_labels()
labels=['Full','-N', '-P']
l = ax.legend(handles[:3], labels, loc=0, handlelength=0.8, borderpad=0.5)


plt.axvline(0.5,color='black',linestyle= '--', linewidth=1.5)
plt.axvline(1.5,color='black',linestyle= '--', linewidth=1.5)

ax.set_xlim(-0.5,2.5)
ax.set_ylim(5,60)
#ax.set_yticklabels([0,0,1,2,3,4,5],size=12,color='k')
ax.yaxis.set_major_formatter(pter)
ax.set_xticklabels(['Zm','Sb','Pv'],size=10,color='k', fontstyle='italic')
ax.set_ylabel('Root length',size=10, color='black')
ax.set_xlabel('')

spList = ['Maize', 'Sorghum', 'Paspalum']
for asp, grp in dfbiomass.groupby('Species'):
    if not asp in spList: continue
    vmax = grp.Length.max() + 2
    vheight = 2
    gap = 3
    plotCap(asp, vmax, vheight, gap)

plt.savefig('rootLength.svg')
plt.savefig('rootLength.png', dpi=300)
