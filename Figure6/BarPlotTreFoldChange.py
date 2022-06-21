import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib import rcParams
import scipy
from scipy import stats
import seaborn as sns

dfall = pd.read_csv('Trehalose_atg12.csv')
print(dfall.head())
dfw22 = dfall.loc[dfall['Strain'] == 'W22']
dfatg12 = dfall.loc[dfall['Strain'] == 'ATG12-2']

def ttest_of_tre(df, strain):
   for condition, grp in df.groupby('Condition'):
      NT = grp.loc[grp['Treatment']=='NT', 'Trehalose_scaled'].values
      ValA = grp.loc[df['Treatment']=='ValA', 'Trehalose_scaled'].values
      tCond, p_valA = scipy.stats.ttest_ind(NT, ValA)
      print(strain, condition, p_valA)
      
ttest_of_tre(dfw22, 'W22')
ttest_of_tre(dfatg12, 'ATG12-2')



rcParams['font.family'] = 'sans-serif'
rcParams['font.size'] = 10


def plotCap(species, vmax, vheight, gap):
    h = vheight + gap
    hh = 2*vheight + gap
    
    xDict = {'3MA':[0, 1],'ValA':[0, 2],'ValA+3MA':[0,3]}
    sigDict = {'3MA':'*','ValA':'****','ValA+3MA':'ns'}
    
    x1, x2 = xDict[species][0],xDict[species][1]
    t1 = sigDict[species]

    plt.plot([x1,x1,x2,x2],[vmax, vmax + vheight, vmax + vheight, vmax], color = 'k', linewidth = 1) 
    ax.text((x1 + x2)* .5, vmax + vheight + 0.1, t1, ha ='center', va='center')
   


#'''
order = ['Control', 'ValA', '3MA', 'ValA+3MA']
fig, axes = plt.subplots(1,2, figsize=(4.3,3.65), sharey=True, tight_layout = True) 

plt.setp(axes[0].spines.values(), linewidth=1.5)
plt.setp(axes[0].spines.values(), linewidth=1.5)
axes[0].spines['top'].set_visible(False)
axes[0].spines['right'].set_visible(False)
axes[0].tick_params(axis='both', which='major', labelsize=10, direction='out', length=3, width=1.5)


plt.setp(axes[1].spines.values(), linewidth=1.5)
plt.setp(axes[1].spines.values(), linewidth=1.5)
axes[1].spines['top'].set_visible(False)
axes[1].spines['right'].set_visible(False)
axes[1].tick_params(axis='both', which='major', labelsize=10, direction='out', length=3, width=1.5)

sns.barplot(x = "Condition", y = 'Trehalose_scaled',  hue = 'Treatment', data=dfw22, hue_order = ['NT','ValA'], edgecolor = 'k', linewidth=1.5, palette = 'Set2', ci="sd",capsize=0.1, alpha = 0.5, ax=axes[0])
#sns.boxplot(x = "Treatment", y = 'DryWeight',  data=df,   fliersize=0, boxprops=dict(alpha = 0.5), width = 0.5, linewidth=1, palette = 'Set2', ax=ax)
sns.stripplot(x = "Condition", y = 'Trehalose_scaled',  hue = 'Treatment',hue_order = ['NT','ValA'], data=dfw22, dodge = True,  jitter = 0.02,  s=3,  palette = 'Set2', ax=axes[0])

sns.barplot(x = "Condition", y = 'Trehalose_scaled',  hue = 'Treatment', data=dfatg12, hue_order = ['NT','ValA'], edgecolor = 'k', linewidth=1.5, palette = 'Set2', ci="sd",capsize=0.1, alpha = 0.5, ax=axes[1])
#sns.boxplot(x = "Treatment", y = 'DryWeight',  data=df,   fliersize=0, boxprops=dict(alpha = 0.5), width = 0.5, linewidth=1, palette = 'Set2', ax=ax)
sns.stripplot(x = "Condition", y = 'Trehalose_scaled',  hue = 'Treatment',hue_order = ['NT','ValA'], data=dfatg12, dodge = True,  jitter = 0.02,  s=3,  palette = 'Set2', ax=axes[1])
axes[0].get_legend().remove()
handles, labels = axes[1].get_legend_handles_labels()
labels=['NT','ValA']
l = axes[1].legend(handles[:2], labels, loc=0, handlelength=0.8, borderpad=0.5)
axes[0].set_ylabel('Trehalose relative abundance')
axes[0].set_xlabel('WT')
axes[1].set_ylabel('')
axes[1].set_xlabel('ATG12-2')
axes[0].set_ylim(0,4.5)
axes[1].set_ylim(0,4.5)
#yticks = [0,1,2,3,4]
#axes[0].set_yticks(yticks)
#axes[0].set_yticks(yticks)
def change_width(ax, new_value) :
    for patch in ax.patches :
        current_width = patch.get_width()
        diff = current_width - new_value

        # we change the bar width
        patch.set_width(new_value)

        # we recenter the bar
        patch.set_x(patch.get_x() + diff * .5)

change_width(axes[0], .3)
change_width(axes[1], .3)

plt.savefig('Trehalose_abundance.png', dpi=500)
plt.savefig('Trehalose_abundance.svg')
'''
vmax = df['DryWeight'].max()
i=1
for condition, grp in dfw22.groupby('Condition'):
	if condition == 'Control':continue
	vvmax = vmax + 0.02 *i
	vheight = 0.03
	gap = 0.3*i
    	#t1 = sigDict[species][condition]
	plotCap(condition, vvmax, vheight, gap)
	i+=3
ax.set_ylabel('Above ground dry biomass',size=10)    
plt.savefig('3MA_dryweight.svg')
plt.savefig("3MA_dryweight.png",dpi = 500)
'''
