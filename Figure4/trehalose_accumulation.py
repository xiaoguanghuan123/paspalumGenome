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
rcParams['text.usetex'] = False
spList = ['Maize','Sorghum','Paspalum']
dfcontent = pd.read_csv('TrehaloseAccumulation.csv' )
met = 'Trehalose'

def plotCap(conditions, vmax, vheight, gap, t1):
    h = vheight + gap
    xDict = {'Full':[-0.2, 0.2],'Ndep':[0.8, 1.2]}
    x1, x2 = xDict[condition][0],xDict[condition][1]
    plt.plot([x1,x1,x2,x2],[vmax, vmax + vheight, vmax + vheight, vmax], color = 'k', linewidth = 1)
    ax.text((x1 + x2)* .5, vmax + vheight + 0.01, t1, ha ='center', va='center')
    
sigDict = {'Maize':{'Full':'***','Ndep':'***'}, 'Sorghum':{'Full':'***','Ndep':'***'}, 'Paspalum':{'Full':'*', 'Ndep':'ns'}}

spDict = {'Maize': 'Z. mays', 'Sorghum':'S. bicolor', 'Paspalum':'P. vaginatum'}
for i, species in enumerate(spList):
    fig, ax = plt.subplots(figsize=(2,3.5), tight_layout = True) 
    plt.setp(ax.spines.values(), linewidth=1.5)
    plt.setp(ax.spines.values(), linewidth=1.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='both', which='major', labelsize=10,direction='out', length=3, width=1.5)
    dfsp = dfcontent.loc[dfcontent['Species'] == species]
    print(dfsp.head())
    Content_Full = dfsp.loc[dfsp['Condition'] == 'Full',['Treatment', met]]
    Content_Ndep = dfsp.loc[dfsp['Condition'] == 'Ndep',['Treatment', met]]
    F_nt = Content_Full.loc[Content_Full['Treatment'] == 'Control', met].values 
    F_valA = Content_Full.loc[Content_Full['Treatment'] == 'ValA', met].values 
    N_nt = Content_Ndep.loc[Content_Ndep['Treatment'] == 'Control', met].values 
    F_valA = Content_Ndep.loc[Content_Ndep['Treatment'] == 'ValA', met].values 
    
    tF, pF = scipy.stats.ttest_ind(F_nt, F_valA)
    tN, pN = scipy.stats.ttest_ind(N_nt, F_valA)
    print(species, 'Full', pF)
    print(species, 'Ndep', pN)
    sns.boxplot(x = "Condition", y = met,  hue = 'Treatment', data=dfsp, fliersize=0, boxprops=dict(alpha = 0.5), width = 0.8, linewidth=1, palette = 'Set2', ax=ax)
    sns.stripplot(x = "Condition", y = met, hue = 'Treatment', data=dfsp, dodge = True,  jitter = 0.02,  s=3,  ax=ax, palette = 'Set2')
    ax.set_ylim(0,0.7)
    ax.set_ylabel('Trehalose relative abundance',size=10)    
    #ax.set_xlabel(spDict[species], fontstyle = 'italic', size =10)
    ax.set_xlabel('')
    ax.get_legend().remove()
    ax.set_xticks([0, 1])
    ax.set_title('$\it{0}$ 21 dap'.format(spDict[species])) #$\it{text you want to show in italics}$'
    ax.set_xticklabels(['Full', '-N'],size=10, color='k')

    for condition, grp in dfsp.groupby('Condition'):
    	vmax = grp[met].max() + 0.02
    	vheight = 0.02
    	gap = 0.3
    	t1 = sigDict[species][condition]
    	plotCap(condition, vmax, vheight, gap, t1)
    
    #plt.show()
    plt.savefig('{0}_{1}_Accumulation.svg'.format(species, met))
    plt.savefig("{0}_{1}_Accumulation.png".format(species, met),dpi = 500)
