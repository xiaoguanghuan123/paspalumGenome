import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib import rcParams

rcParams['font.family'] = 'sans-serif'
rcParams['font.size'] = 10


def pcttck(x,pos):
    return "{0}%".format(int(x))
pter = FuncFormatter(pcttck)

nDec, pDec = [], []
nStd, pStd = [], []
spList = ['Maize','Sorghum','Paspalum']


dfcontent = pd.read_csv('DryWeightContent.csv' )
dfcontent = dfcontent.loc[dfcontent['Treatment'] == 'NT']
dfcontent.loc[dfcontent['Condition']=='Full', 'Condition'] = 'Control'
dfcontent.loc[dfcontent['Condition']=='Ndep', 'Condition'] = '-N'
dfcontent.loc[dfcontent['Condition']=='Pdep', 'Condition'] = '-P'
dfN = dfcontent.drop('P', axis = 1)
dfP = dfcontent.drop('N', axis = 1)
print(list(dfN))
dfreduction = pd.DataFrame()
 

def calculate_mean_std(df, condition, meas):
   ctr = df.loc[df['Condition'] == 'Control', meas]
   dep = df.loc[df['Condition'] == condition, meas]
   redList = list(map(lambda x, y: 100*(x - y)/x, ctr, dep))
   redMean = np.mean(redList)
   redStd = np.std(redList)
   return redMean, redStd
   
for species in spList:
    dfspP = dfP.loc[dfP['Species'] == species]
    dfspN = dfN.loc[dfN['Species'] == species]
    Nred_mean, Nred_std = calculate_mean_std(dfspN, '-N', 'N')
    Pred_mean, Pred_std = calculate_mean_std(dfspP, '-P', 'P')
    nDec.append(Nred_mean)
    nStd.append(Nred_std)
    pDec.append(Pred_mean)
    pStd.append(Pred_std)
print(nDec, nStd)
print(pDec, pStd)
#'''          
ind = [0,1,2]
width = 0.8
fig, axes = plt.subplots(1,2, figsize=(5.5,3), tight_layout = True)   
plt.setp(axes[0].spines.values(), linewidth=1.5)
plt.setp(axes[1].spines.values(), linewidth=1.5)
#axes[0].spines['top'].set_visible(False)
axes[0].spines['right'].set_visible(False)
axes[0].tick_params(axis='both', which='major', labelsize=10,direction='out', length=4, width=1.5)
#axes[1].spines['top'].set_visible(False)
axes[1].spines['right'].set_visible(False)
axes[1].tick_params(axis='both', which='major', labelsize=10,direction='out', length=4, width=1.5)


axes[0].bar(ind, nDec, yerr = nStd, capsize=5, width = 0.5, edgecolor = 'k', color = 'orange')
axes[0].invert_yaxis()
axes[1].bar(ind, pDec, yerr = nStd, capsize=5, width = 0.5, edgecolor = 'k', color = 'green')
axes[1].invert_yaxis()
axes[0].yaxis.set_major_formatter(pter)
axes[1].yaxis.set_major_formatter(pter)

axes[0].set_xticks(ind)
axes[0].set_xticklabels(['Zm','Sb','Pv'],size=10,color='k', fontstyle='italic' )
axes[1].set_xticks(ind)
axes[1].set_xticklabels(['Zm','Sb','Pv'],size=10,color='k', fontstyle='italic')

axes[0].set_ylabel('Reduction in N')
axes[1].set_ylabel('Reduction in P')

plt.savefig('PercentageOfReduction.svg')
plt.savefig("PercentageOfReduction.png",dpi = 500)
#'''
