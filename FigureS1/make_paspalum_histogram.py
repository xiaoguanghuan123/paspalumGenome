import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.family'] = 'sans-serif'
rcParams['font.size'] = 10

pseudo_covs = []
scaff_covs = []
fh = open("paspalumcoverage.csv")
for x in fh:
    y = x.strip().split(',')
    if len(y) < 5: break
    if float(y[4]) > 120: continue
    if "scaffold" in y[0]:
        scaff_covs.append(float(y[4]))
    else:
        pseudo_covs.append(float(y[4]))
figure = plt.figure(figsize=(4,3.8), tight_layout=True)
ax = figure.add_subplot(1,1,1)

plt.setp(ax.spines.values(), linewidth=1.5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', which='major', labelsize=10, direction='out', length=3, width=1.5)

ax.hist([pseudo_covs],bins=120,stacked=True, edgecolor ='k', linewidth=0.5, label=["Pseudomolecules"])
ax.legend()
ax.set_xlabel("Average sequencing depth")
ax.set_ylabel("Frequency (among 10 kb windows)")
plt.savefig('SeqCov.svg', dpi=500)
