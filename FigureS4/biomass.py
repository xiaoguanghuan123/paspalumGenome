import pandas as pd
import seaborn as sns
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib import rcParams
import numpy as np

def pcttck(x,pos):
    return "{0}g".format(int(x) )
pter = FuncFormatter(pcttck)

dfbiomass = pd.read_csv('timepoint_dryweight.csv')
#dfbiomass.loc[dfbiomass['Condition']=='Full', 'Condition'] = 'Control'

group = dfbiomass.groupby('DAP')
for g in group:
    #print(g[1])
    print(g[0])
    g_sorg = g[1].loc[g[1]['Species'] == 'Sorghum']
    Full = g_sorg.loc[g_sorg['Condition']=='Full', 'Biomass']
    N_dep = g_sorg.loc[g_sorg['Condition']=='N', 'Biomass']
    tFN = ttest_ind(Full,N_dep)
    print('sorghum','N_dep:', tFN[1])
    
    g_pas = g[1].loc[g[1]['Species'] == 'Paspalum']
    Full = g_pas.loc[g_pas['Condition']=='Full', 'Biomass']
    N_dep = g_pas.loc[g_pas['Condition']=='N', 'Biomass']
    tFN = ttest_ind(Full,N_dep)
    print('paspalum','N_dep:', tFN[1])

    


rcParams['font.family'] = 'sans-serif'
rcParams['font.size'] = 10





for species, grp in dfbiomass.groupby('Species'):
	grp = grp.loc[grp['DAP'] != 9]
	fig, ax = plt.subplots(figsize=(2.8,3), tight_layout = True)
	plt.setp(ax.spines.values(), linewidth=1.5)
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.tick_params(axis='both', which='major', labelsize=10,direction='out', length=4, width=1.5)
	sns.barplot( x="DAP", hue="Condition", y="Biomass", data=grp,ax=ax)
	sns.stripplot(x="DAP", y="Biomass", hue="Condition", data=grp, dodge=True, linewidth=1, s=2, ax=ax)
	#ax.get_legend().remove()
	#handles, labels = ax.get_legend_handles_labels()
	#l = ax.legend(handles[:3], labels[:3])
	handles, labels = ax.get_legend_handles_labels()
	labels=['Full','-N', '-P']
	l = ax.legend(handles[:3], labels, loc=0, handlelength=0.8, borderpad=0.5)
	ax.set_title(species)
	plt.savefig('Biomass_{0}.png'.format(species), dpi=300)
	plt.show()




