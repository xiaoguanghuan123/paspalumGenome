import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import matplotlib as mpl
from matplotlib.ticker import FuncFormatter
from matplotlib import rcParams
from scipy.stats import ttest_ind



def lentck(x, pos):
    return "{0}cm".format(int(x))
lenter = FuncFormatter(lentck)

def fdstck(x, pos):
    return "{0}dap".format(int(x))
fdster = FuncFormatter(fdstck)

def masstck(x, pos):
    return "{0}g".format(int(x))
masster = FuncFormatter(masstck)

rcParams['font.family']=['sans-serif']
rcParams['font.size']=10



df = pd.read_csv('AdultPlantPhenotypeForStat.csv')


        
fig,axes = plt.subplots(1,4, figsize=(8,2.6), tight_layout=True)
for a, b in enumerate(list(df)[2:]):
    F = df.loc[df['Condition'] == 'Full',b]
    Fv = df.loc[df['Condition'] == 'FV',b]
    tFN = ttest_ind(F,Fv)
    print(b,'FV:', tFN[1])
    ax = axes[a]
    plt.setp(ax.spines.values(), linewidth=1.5)
    plt.setp(ax.spines.values(), linewidth=1.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='both', which='major', labelsize=10, direction='out', length=3, width=1.5)
    title = b.split('_')[0]
    sns.boxplot(x = 'Condition', y = b, data = df, fliersize=0, boxprops=dict(alpha = 0.5), width = 0.5, linewidth=1, palette = 'Set2',ax=ax)
    sns.stripplot(x='Condition',y=b,data=df, jitter=True, dodge = True,  s=5,  ax=ax, palette = 'Set2')
    ax.set_xticklabels(['Control','ValA'],size=10)
    ax.set_ylabel(title,size=10, color='black')
    ax.set_xlabel('')
    if a < 2 : 
    	ax.yaxis.set_major_formatter(lenter)
    	t1 = '*'
    if a == 2 : 
    	ax.yaxis.set_major_formatter(fdster)
    	t1 = '***'
    if a == 3 : 
    	ax.yaxis.set_major_formatter(masster)
    	t1 = 'ns'
    	 
    x1, x2 = 0, 1
    vmax = df[b].max() * 1.05
    vheight = df[b].max() * 0.02
    ax.plot([x1,x1,x2,x2],[vmax, vmax + vheight, vmax + vheight, vmax], color = 'k', linewidth = 1)
    ax.text((x1 + x2)* .5, (vmax + vheight)*1.01, t1, ha ='center', va='center')
    plt.savefig('{0}.svg'.format(title))
