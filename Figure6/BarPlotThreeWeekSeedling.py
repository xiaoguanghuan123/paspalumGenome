import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib import rcParams
import scipy
from scipy import stats
import seaborn as sns

df = pd.read_csv('ThreeWeekSeedling_3MA.csv')
df.head()

nt = df.loc[df['Treatment']=='Control', 'DryWeight'].values
MA = df.loc[df['Treatment']=='3MA', 'DryWeight'].values
valA = df.loc[df['Treatment']=='ValA', 'DryWeight'].values
double = df.loc[df['Treatment']=='ValA+3MA', 'DryWeight'].values

tF, p_3MA = scipy.stats.ttest_ind(nt, MA)
tF, p_valA = scipy.stats.ttest_ind(nt, valA)
tF, p_double = scipy.stats.ttest_ind(nt, double)
tF, p_valA_double = scipy.stats.ttest_ind(MA, double)

print('3MA: ', p_3MA)
print('valA: ', p_valA)
print('double: ',p_double)
print('MA vs double: ', p_valA_double)

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
   


def pcttck(x,pos):
    return "{0}g".format(str(float(x))[:3])
pter = FuncFormatter(pcttck)

order = ['Control', 'ValA', '3MA', 'ValA+3MA']
fig, ax = plt.subplots(figsize=(3.5,3.7), tight_layout = True) 
plt.setp(ax.spines.values(), linewidth=1.5)
plt.setp(ax.spines.values(), linewidth=1.5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', which='major', labelsize=10, direction='out', length=3, width=1.5)
sns.barplot(x = "Treatment", y = 'DryWeight',  data=df, order = order, edgecolor = 'k', linewidth=1, palette = 'Set2', ci="sd",capsize=0.1,ax=ax)
#sns.boxplot(x = "Treatment", y = 'DryWeight',  data=df,   fliersize=0, boxprops=dict(alpha = 0.5), width = 0.5, linewidth=1, palette = 'Set2', ax=ax)
sns.stripplot(x = "Treatment", y = 'DryWeight',  data=df,order = order, dodge = True,  jitter = 0.02,  s=3,  ax=ax, palette = 'Set2')
#ax.set_xlabel('Z. mays', fontstyle = 'italic')
ax.set_xlabel('')
ax.set_xticks([0,1,2,3])
#ax.set_title('WT 21 dap')
ax.set_xticklabels(['Control', 'ValA','3-MA', 'ValA+3-MA'],rotation = 20, ha = 'right')
ax.set_ylim(0,1)

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

vmax = df['DryWeight'].max()
i=1
for condition, grp in df.groupby('Treatment'):
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
