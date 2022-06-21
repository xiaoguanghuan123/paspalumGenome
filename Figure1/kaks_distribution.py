import os 
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


dfka = pd.read_csv('SevenSpeciesKA_V3.csv', sep = '\t', index_col = 0)
dfks = pd.read_csv('SevenSpeciesKS_V3.csv',sep='\t', index_col = 0)
dfkaks = pd.read_csv('SevenSpeciesKAKS_V3.csv',sep='\t', index_col = 0)
dfka = dfka.loc[dfkaks.index]
dfks = dfks.loc[dfkaks.index]

dfSynMerge = pd.read_csv('Syntenic_GENE_SET.csv', index_col = 0)

plotOrder = ['Zmays','Sobic','Pavag','Seita', 'Ortho', 'Bradi', 'Orysa']





import matplotlib as mpl
import matplotlib.font_manager as font_manager
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.size'] = 12

sns.set(style = 'darkgrid')
mpl.rcParams['font.family']=['sans-serif']
fig, axes = plt.subplots(1,2,figsize=(12,6), tight_layout = True)

dfks[plotOrder].plot.kde(xlim=(0,1.2),linewidth=2,ax=axes[0])
axes[0].set_xlabel('Synonymous nucleotide substituion rate')
font = font_manager.FontProperties(family='sans-serif', style='italic', size=12)                                  
ksLegend = axes[0].legend(labels = ['Z. mays','S. bicolor','P. vaginatume', 'S. italica', 'O. thomaeum', 'B. distathyon', 'O. sativa'], prop = font)
for line in ksLegend.get_lines():
    line.set_linewidth(4.0)

dfkaks = dfkaks[plotOrder]
dfkaksmelt = dfkaks.melt(var_name='Species', value_name='KA/KS')

sns.boxplot(x = "Species", hue = "Species", y = "KA/KS",  dodge=False, data=dfkaksmelt,  ax=axes[1],fliersize=2)
sns.pointplot(x='Species', y="KA/KS", data=dfkaksmelt.groupby('Species', as_index=True).median().loc[plotOrder].reset_index(), color = 'k' ,ax=axes[1])
axes[1].get_legend().remove()
axes[1].set_xticklabels(['Z. mays','S. bicolor','P. vaginatume', 'S. italica', 'O. thomaeum', 'B. distathyon', 'O. sativa'],size=12,color='k', fontstyle='italic', rotation=20, ha = 'right')
axes[1].set_title('Ka/Ks value in all orghologous genes',size=12)

axes[1].set_xlabel('')
axes[1].set_ylabel('Ka/Ks',size=12)


fig.savefig('KAKSComparison.svg')
fig.savefig('KAKSComparison.png',dpi=150)
