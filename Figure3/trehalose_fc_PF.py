import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

fcDir = 'FoldChange_ttest_GS/'
fcList = [afile for afile in os.listdir(fcDir) if afile.startswith('Combined')]
fcNList = [os.path.join(fcDir, afile) for afile in fcList if 'Ndep' in afile]
fcPList = [os.path.join(fcDir, afile) for afile in fcList if 'Pdep' in afile]

dffc_P = pd.DataFrame()
treDict = {}
for i, afile in enumerate(fcPList):
    sp = afile.split('/')[1].split('_')[2]
    dfPfc = pd.read_csv(afile, sep = '\t', index_col = 0)
    dfPfc = dfPfc.iloc[:,[1,2,4]]
    treDict[sp]= {}
    treDict[sp]['log2fc'] = dfPfc.loc['Trehalose']['log2(FC)']
    treDict[sp]['pval'] = dfPfc.loc['Trehalose']['p.value']
    treDict[sp]['condition'] = 'Pdep' 
  
    #print(dfNfc.head())
    foldChange = []
    for x in dfPfc.index:
        if dfPfc.loc[x]['p.value']>0.05: foldChange.append(0)
        elif abs(dfPfc.loc[x]['log2(FC)']) <1 : foldChange.append(0)
        else: foldChange.append(dfPfc.loc[x]['log2(FC)'])
    dfPfc[sp] = foldChange
    if i == 0: 
        dffc_P = pd.DataFrame(dfPfc[sp])
    else:
        dfPfc = dfPfc.loc[dffc_P.index]
        dffc_P[sp] = dfPfc[sp]
        
        
for x in dffc_P.index:
    i=0
    for v in dffc_P.loc[x].values:
        if abs(v) <=1: i += 1
    if i == 3: dffc_P = dffc_P.drop(x, axis = 0)

nsList = []
for x in dffc_P.index:
    vCount = dict(dffc_P.loc[x].value_counts())
    if 0 in vCount: nsList.append(vCount[0])
    else: nsList.append(0)

dffc_P['nsCount'] = nsList
dffc_P.sort_values(by = 'nsCount', ascending = True, inplace=True)
dffc_P = dffc_P.drop('nsCount',axis =1)
dftre_P = pd.DataFrame.from_dict(treDict,orient = 'index')
dftre_P['species'] = ['Pv','Zm', 'Sb']

fig, ax = plt.subplots(figsize=(1.7,2), tight_layout = True) 
plt.setp(ax.spines.values(), linewidth=1.5)
plt.setp(ax.spines.values(), linewidth=1.5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', which='major', labelsize=10, direction='out', length=3, width=1.5)
spDict={'Sb':{'x':3, 'c':'grey'},'Pv':{'x':1, 'c':'purple'},'Zm':{'x':2, 'c':'grey'} }
for sp, df in dftre_P.groupby('species'):
    ax.bar(x=spDict[sp]['x'], height = df['log2fc'], color = spDict[sp]['c'], width = 0.7, edgecolor = 'k')

#sns.barplot(x= 'condition', y ='log2fc', hue = 'species', edgecolor = 'k', 
#            linewidth=1, palette = 'tab20_r', data = dftre_P, ax = ax) 

xticks = [1,2,3]
def change_width(ax, new_value) :
    for patch in ax.patches :
        current_width = patch.get_width()
        diff = current_width - new_value
        current_height = patch.get_height()
        # we change the bar width
        patch.set_width(new_value)
        
        # we recenter the bar
        newcenter = patch.get_x() + diff * .5
        newxtick = patch.get_x() + diff*2
        print(patch.get_x())
        xticks.append(newxtick)
        patch.set_x(newcenter)
        
#change_width(ax, .2)

ax.set_ylim(-0.5,3)        
#ax.get_legend().remove()
ax.set_xlabel('')
ax.set_xticks(xticks)
ax.set_xticklabels(['Pv','Zm', 'Sb'], fontstyle = 'italic')    
ax.set_ylabel('Trehalose log$r_{2}$(-P/Full)',size=10)    #$\it{text you want to show in italics}$
#ax.set_title('$\it{0}$'.format(geneName)) 
plt.savefig('Trehalose_fc_Pdep.svg')
plt.savefig('Trehalose_fc_Pdep.png', dpi=300)
