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

spList = ['Maize','Sorghum','Paspalum']
spDict = {'Maize': 'Z. mays', 'Sorghum':'S. bicolor', 'Paspalum':'P. vaginatum'}
dfbiomass = pd.read_csv('DryWeight_ValA.csv')
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
    
sigDict = {'Maize':{'Full':'**','-N':'*', '-P':'**'}, 'Sorghum':{'Full':'ns','-N':'*', '-P':'ns'}, 'Paspalum':{'Full':'ns', '-N':'ns', '-P':'ns'}}
yDict = {'Maize':1, 'Sorghum':0.3, 'Paspalum':0.5}
for i, species in enumerate(spList):
    fig, ax = plt.subplots(figsize=(2,2.7), tight_layout = True) 
    plt.setp(ax.spines.values(), linewidth=1.5)
    plt.setp(ax.spines.values(), linewidth=1.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='both', which='major', labelsize=10,direction='out', length=3, width=1.5)
    dfsp = dfbiomass.loc[dfbiomass['Species'] == species]
    print(dfsp.head())
    Content_Full = dfsp.loc[dfsp['Condition'] == 'Full',['Treatment', met]]
    Content_Ndep = dfsp.loc[dfsp['Condition'] == '-N',['Treatment', met]]
    Content_Pdep = dfsp.loc[dfsp['Condition'] == '-P',['Treatment', met]]
    
    F_nt = Content_Full.loc[Content_Full['Treatment'] == 'Control', met].values 
    F_valA = Content_Full.loc[Content_Full['Treatment'] == 'ValA', met].values 
    N_nt = Content_Ndep.loc[Content_Ndep['Treatment'] == 'Control', met].values 
    N_valA = Content_Ndep.loc[Content_Ndep['Treatment'] == 'ValA', met].values 
    P_nt = Content_Ndep.loc[Content_Ndep['Treatment'] == 'Control', met].values 
    P_valA = Content_Ndep.loc[Content_Ndep['Treatment'] == 'ValA', met].values
    tF, pF = scipy.stats.ttest_ind(F_nt, F_valA)
    tN, pN = scipy.stats.ttest_ind(N_nt, N_valA)
    tP, pP = scipy.stats.ttest_ind(P_nt, P_valA)
    print(species, 'Full', pF)
    print(species, '-N', pN)
    print(species, '-P', pF)
    dfsp = dfsp.loc[dfsp['Condition']!='-P']
    sns.boxplot(x = "Condition", y = met,  hue = 'Treatment', data=dfsp, fliersize=0, boxprops=dict(alpha = 0.5), width = 0.8, linewidth=1,palette = 'Set2',  ax=ax)
    sns.stripplot(x = "Condition", y = met, hue = 'Treatment', data=dfsp, dodge = True,  jitter = 0.02,  s=3, palette = 'Set2',  ax=ax)
    ylim = yDict[species]
    ax.set_ylim(0,ylim)
    ax.set_ylabel('Above ground dry biomass',size=10)    
    ax.get_legend().remove()
    
    ax.set_xticklabels(['Full', '-N'],size=10, color='k')
    ax.yaxis.set_major_formatter(pter)
    ax.set_xlabel(spDict[species], fontstyle = 'italic', size =10)
    for condition, grp in dfsp.groupby('Condition'):
    	vmax = grp[met].max() + 0.02
    	vheight = 0.02
    	gap = 0.3
    	t1 = sigDict[species][condition]
    	plotCap(condition, vmax, vheight, gap, t1)
    
    
    plt.savefig('{0}_{1}_dryweight.svg'.format(species, met))
    plt.savefig("{0}_{1}_dryweight.png".format(species, met),dpi = 500)
    
    
'''    
dfAbove = dfbiomass.drop(labels= ['Root','Shoot_to_root'],axis=1)
dfRoot = dfbiomass.drop(labels=['AboveGround','Shoot_to_root'],axis=1)
dfSR = dfbiomass.drop(labels=['Root','AboveGround'],axis=1)
GroupAbove = dfAbove.groupby('Species')
GroupRoot = dfRoot.groupby('Species')
GroupSR = dfSR.groupby('Species')

for i,g in enumerate(GroupAbove):
    ax = axes[0,i]
    Species = list(set(g[1]['Species']))[0]
    print(Species)
    sns.boxplot(x="Condition",  y="AboveGround", hue="Treatment", data=g[1], dodge=True, fliersize=0,linewidth=3,ax=ax)
    sns.stripplot(x="Condition", y="AboveGround", hue="Treatment", data=g[1], dodge=True,linewidth=3, ax=ax)
    
    handles, labels = ax.get_legend_handles_labels()
    l = ax.legend(handles[:2], labels[:2])
    
    ax.set_xlim(-0.5,2.5)
    #ax.set_ylim(0,1) #for dry weight
    ax.set_ylim(0,1) #for fresh weight
    ax.axvline(0.5,color='black',linestyle= '--', linewidth=1.5)
    ax.axvline(1.5,color='black',linestyle= '--', linewidth=1.5)
    ax.set_title(Species)
    dfGroupAbove = g[1].set_index(['Treatment','Condition'])
    Full_NO_ValA = dfGroupAbove.loc['NO_ValA'].loc['Full'].AboveGround.values
    N_dep_NO_ValA = dfGroupAbove.loc['NO_ValA'].loc['-N'].AboveGround.values
    P_dep_NO_ValA = dfGroupAbove.loc['NO_ValA'].loc['-P'].AboveGround.values
    Full_ValA = dfGroupAbove.loc['ValA'].loc['Full'].AboveGround.values
    N_dep_ValA = dfGroupAbove.loc['ValA'].loc['-N'].AboveGround.values
    P_dep_ValA = dfGroupAbove.loc['ValA'].loc['-P'].AboveGround.values
    tFULL = scipy.stats.ttest_ind(Full_NO_ValA,Full_ValA)
    tNdep = scipy.stats.ttest_ind(N_dep_NO_ValA,N_dep_ValA)
    tPdep = scipy.stats.ttest_ind(P_dep_NO_ValA,P_dep_ValA)
    #tFP_ValA = ttest_ind(Full_ValA,P_dep_ValA)
    print(Species,'Full:', tFULL[1])
    print(Species,'N_dep:', tNdep[1])
    print(Species,'P_dep:', tPdep[1])
    
for i,g in enumerate(GroupSR):
    ax = axes[1,i]
    Species = list(set(g[1]['Species']))[0]
    print(Species)
    sns.boxplot(x="Condition",  y="Shoot_to_root", hue="Treatment", data=g[1], dodge=True, fliersize=0,linewidth=3,ax=ax)
    sns.stripplot(x="Condition", y="Shoot_to_root", hue="Treatment", data=g[1], dodge=True,linewidth=3, ax=ax)
    
    handles, labels = ax.get_legend_handles_labels()
    l = ax.legend(handles[:2], labels[:2])
    
    ax.set_xlim(-0.5,2.5)
    ax.set_ylim(0,15)
    ax.axvline(0.5,color='black',linestyle= '--', linewidth=1.5)
    ax.axvline(1.5,color='black',linestyle= '--', linewidth=1.5)
    
    dfGroupSR = g[1].set_index(['Treatment','Condition'])
    Full_NO_ValA = dfGroupSR.loc['NO_ValA'].loc['Full'].Shoot_to_root.values
    N_dep_NO_ValA = dfGroupSR.loc['NO_ValA'].loc['-N'].Shoot_to_root.values
    P_dep_NO_ValA = dfGroupSR.loc['NO_ValA'].loc['-P'].Shoot_to_root.values
    Full_ValA = dfGroupSR.loc['ValA'].loc['Full'].Shoot_to_root.values
    N_dep_ValA = dfGroupSR.loc['ValA'].loc['-N'].Shoot_to_root.values
    P_dep_ValA = dfGroupSR.loc['ValA'].loc['-P'].Shoot_to_root.values
    tFULL = scipy.stats.ttest_ind(Full_NO_ValA,Full_ValA)
    tNdep = scipy.stats.ttest_ind(N_dep_NO_ValA,N_dep_ValA)
    tPdep = scipy.stats.ttest_ind(P_dep_NO_ValA,P_dep_ValA)
    #tFP_ValA = ttest_ind(Full_ValA,P_dep_ValA)
    print(Species,'Full:', tFULL[1])
    print(Species,'N_dep:', tNdep[1])
    print(Species,'P_dep:', tPdep[1])
    
    
'''


