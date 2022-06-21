import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib import rcParams
import scipy
from scipy import stats
import seaborn as sns
rcParams['font.family'] = 'sans-serif'
rcParams['font.size'] = 10

spList = ['W22','ATG12-2']
spDict = {'Maize': 'Z. mays', 'Sorghum':'S. bicolor', 'Paspalum':'P. vaginatum'}
dfbiomass = pd.read_csv('W22_ATG12-2_DryWeight_ValA.csv')
met = 'AboveGround'
def plotCap(conditions, vmax, vheight, gap, t1):
    h = vheight + gap
    xDict = {'Full':[-0.2, 0.2],'-N':[0.8, 1.2], '-P':[1.8, 2.2]}
    x1, x2 = xDict[condition][0],xDict[condition][1]
    plt.plot([x1,x1,x2,x2],[vmax, vmax + vheight, vmax + vheight, vmax], color = 'k', linewidth = 1)
    ax.text((x1 + x2)* .5, vmax + vheight + 0.01, t1, ha ='center', va='center')
def pcttck(x,pos):
    return "{0}g".format(str(x)[:3] )
pter = FuncFormatter(pcttck)
    
sigDict = {'W22':{'Full':'**','-N':'*'}, 'ATG12-2':{'Full':'ns','-N':'*', '-P':'ns'}}
yDict = {'Maize':1, 'Sorghum':0.3, 'Paspalum':0.5}


fig, axes = plt.subplots(1,2, figsize=(4.5,3.5), tight_layout = True, sharey=True) 
for i, species in enumerate(spList):
    ax = axes[i]
    plt.setp(ax.spines.values(), linewidth=1.5)
    plt.setp(ax.spines.values(), linewidth=1.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='both', which='major', labelsize=10,direction='out', length=3, width=1.5)
    dfsp = dfbiomass.loc[dfbiomass['Species'] == species]
    print(dfsp.head())
    Content_Full = dfsp.loc[dfsp['Condition'] == 'Full',['Treatment', 'AboveGround']]
    Content_Ndep = dfsp.loc[dfsp['Condition'] == '-N', ['Treatment','AboveGround']]
   
    
    F_nt = Content_Full.loc[Content_Full['Treatment'] == 'NT', 'AboveGround'].values 
    F_valA = Content_Full.loc[Content_Full['Treatment'] == 'ValA', 'AboveGround'].values 
    
    
    N_nt = Content_Ndep.loc[Content_Ndep['Treatment'] == 'NT', 'AboveGround'].values 
    N_valA = Content_Ndep.loc[Content_Ndep['Treatment'] == 'ValA', 'AboveGround'].values 
    
    tF, pF = scipy.stats.ttest_ind(F_nt, F_valA)
    tN, pN = scipy.stats.ttest_ind(N_nt, N_valA)
    
    print(species, 'Full', pF)
    print(species, '-N', pN)

    sns.boxplot(x = "Condition", y = met,  hue = 'Treatment', data=dfsp, fliersize=0, boxprops=dict(alpha = 0.5), width = 0.7, linewidth=1,palette = 'Set2',  ax=ax)
    sns.stripplot(x = "Condition", y = met, hue = 'Treatment', data=dfsp, dodge = True,  jitter = 0.02,  s=3, palette = 'Set2',  ax=ax)
    #ylim = yDict[species]
    #ax.set_ylim(0,ylim)
    axes[0].set_ylabel('Above ground dry biomass',size=10)    
    
    
    ax.set_xticklabels(['Full', '-N'],size=10, color='k')
    ax.yaxis.set_major_formatter(pter)
    ax.set_title(species, fontstyle = 'italic', size =10)
    ax.set_ylim(0, 2.2)
    ax.set_xlabel('')
axes[0].get_legend().remove()
axes[1].set_ylabel('')
handles, labels = axes[1].get_legend_handles_labels()
labels=['NT','ValA']
l = axes[1].legend(handles[:2], labels, loc=0, handlelength=0.8, borderpad=0.5)


plt.savefig('W22_ATG12_dryweight.svg')
plt.savefig("W22_ATG12_dryweight.png",dpi = 500)    
#    plt.show()


