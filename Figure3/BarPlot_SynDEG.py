from matplotlib import pyplot as plt
import matplotlib as mpl
import seaborn as sns
import pandas as pd
mpl.rcParams.update(mpl.rcParamsDefault)

mpl.rcParams['axes.linewidth'] = 2
mpl.rcParams['font.family']=['sans-serif']

DeDict={'Zmays':{'WholeGenome':{'NonSyn':18604,'Syn':27737},
                 'NF_Moderate':{'NonSyn':870,'Syn':2187},
                 'PF_Moderate':{'NonSyn':156,'Syn':435}},
       'Sbicolor':{'WholeGenome':{'NonSyn':29027,'Syn':24596},
                   'NF_Moderate':{'NonSyn':1084,'Syn':2060},
                   'PF_Moderate':{'NonSyn':915,'Syn':1403}},
       'Pvaginatum':{'WholeGenome':{'NonSyn':14527,'Syn':17287}, #41,771 annotated genes on chromosome
                   'NF_Moderate':{'NonSyn':1199,'Syn':1524},
                   'PF_Moderate':{'NonSyn':710,'Syn':988}}}

dfZm = pd.DataFrame.from_dict(DeDict['Zmays'],orient='index')
dfSb = pd.DataFrame.from_dict(DeDict['Sbicolor'],orient='index')
dfPv = pd.DataFrame.from_dict(DeDict['Pvaginatum'],orient='index')
df = pd.concat([dfZm,dfSb,dfPv])
df['Species'] = ['Zmays','Zmays','Zmays','Sbicolor','Sbicolor','Sbicolor','Pvaginatum','Pvaginatum','Pvaginatum']
df.reset_index(inplace=True)

dfZm = pd.DataFrame.from_dict(DeDict['Zmays'],orient='index')
dfSb = pd.DataFrame.from_dict(DeDict['Sbicolor'],orient='index')
dfPv = pd.DataFrame.from_dict(DeDict['Pvaginatum'],orient='index')

dfstack_NF = df.set_index('index').loc[['NF_Moderate']]
dfstack_PF = df.set_index('index').loc[['PF_Moderate']]
sns.set_style("white")
sns.set_context({"figure.figsize": (5, 6)})
dfstack_NF['total'] = dfstack_NF.NonSyn + dfstack_NF.Syn
dfstack_PF['total'] = dfstack_PF.NonSyn + dfstack_PF.Syn

from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.size'] = 10

plotOrder = ['Pvaginatum', 'Sbicolor','Zmays']
fig, (ax1,ax2)= plt.subplots(1,2,sharey=True, figsize=(3,3), tight_layout=True)
plt.setp(ax1.spines.values(), linewidth=1.5)
plt.setp(ax2.spines.values(), linewidth=1.5)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax1.tick_params(axis='both', which='major', labelsize=10, direction='out', length=3.5, width=1.5)
ax2.tick_params(axis='both', which='major', labelsize=10, direction='out', length=3.5, width=1.5)

sns.barplot(x = dfstack_NF.Species, y = dfstack_NF.total, color = "pink", edgecolor ='k',  order = plotOrder, linewidth=1, ax=ax1, )
sns.barplot(x = dfstack_NF.Species, y = dfstack_NF.Syn, color = "lightblue",edgecolor ='k', order = plotOrder,linewidth=1, ax=ax1)
sns.barplot(x = dfstack_PF.Species, y = dfstack_PF.total, color = "pink",edgecolor ='k', order = plotOrder,linewidth=1, ax=ax2,label = 'Nonsyntenic')
sns.barplot(x = dfstack_PF.Species, y = dfstack_PF.Syn, color = "lightblue",edgecolor ='k', order = plotOrder,linewidth=1,ax=ax2, label = 'Syntenic')
#ax1.set_xticklabels=(['Zm','Sb','Pv'])
ax2.set_ylabel('')
ax1.set_xlabel('')
ax2.set_xlabel('')
ax1.set_ylabel('Number of DEGs')
ax1.set_xticklabels(['Pv','Zm','Sb'],size=10,color='k', fontstyle='italic')
ax2.set_xticklabels(['Pv','Zm','Sb'],size=10,color='k', fontstyle='italic')
ax1.set_title('DEG(-N)')
ax2.set_title('DEG(-P)')

#ax2.legend()
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
        #print(patch.get_x())
        #xticks.append(newxtick)
        patch.set_x(newcenter)
change_width(ax1,.6)
change_width(ax2,.6)
plt.tight_layout()
plt.savefig('DEG_SYN_CLASS.svg')
