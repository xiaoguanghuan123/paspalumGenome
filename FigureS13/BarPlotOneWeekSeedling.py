import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib import rcParams
import scipy
from scipy import stats
import seaborn as sns

df = pd.read_csv('oneWeekSeedling.csv')
df.head()

nt = df.loc[df['Treatment']=='Control', 'DryWeight'].values
valA = df.loc[df['Treatment']=='ValA', 'DryWeight'].values
tF, pF = scipy.stats.ttest_ind(nt, valA)
print(pF)

rcParams['font.family'] = 'sans-serif'
rcParams['font.size'] = 10

def pcttck(x,pos):
    return "{0}g".format(str(float(x))[:3])
pter = FuncFormatter(pcttck)

fig, ax = plt.subplots(figsize=(2.5,3.4), tight_layout = True) 
plt.setp(ax.spines.values(), linewidth=1.5)
plt.setp(ax.spines.values(), linewidth=1.5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', which='major', labelsize=10, direction='out', length=3, width=1.5)
sns.barplot(x = "Treatment", y = 'DryWeight',  data=df,  edgecolor = 'k', linewidth=1, palette = 'Set2', ci="sd",capsize=0.2,ax=ax)
#sns.boxplot(x = "Treatment", y = 'DryWeight',  data=df,   fliersize=0, boxprops=dict(alpha = 0.5), width = 0.5, linewidth=1, palette = 'Set2', ax=ax)
sns.stripplot(x = "Treatment", y = 'DryWeight',  data=df, dodge = True,  jitter = 0.02,  s=3,  ax=ax, palette = 'Set2')
#ax.set_xlabel('Z. mays', fontstyle = 'italic')
ax.set_xlabel('')
ax.set_xticks([0,1])
ax.set_xticklabels(['NT', 'ValA'])
ax.set_ylim(0,0.3)

ax.yaxis.set_major_formatter(pter)

def change_width(ax, new_value) :
    for patch in ax.patches :
        current_width = patch.get_width()
        diff = current_width - new_value

        # we change the bar width
        patch.set_width(new_value)

        # we recenter the bar
        patch.set_x(patch.get_x() + diff * .5)

change_width(ax, .5)

x1, x2 = 0,1
vmax = df['DryWeight'].max() + 0.02
vheight = 0.02
t1 = 'ns'
plt.plot([x1,x1,x2,x2],[vmax, vmax + vheight, vmax + vheight, vmax], color = 'k', linewidth = 1)
ax.text((x1 + x2)* .5, vmax + vheight + 0.02, t1, ha ='center', va='center')

ax.set_ylabel('Above ground dry biomass',size=10)   
ax.set_xlabel('') 
plt.savefig('OneWeekOldSeedling_dryweight.svg')
plt.savefig("OneWeekOldSeedling_dryweight.png",dpi = 500)
