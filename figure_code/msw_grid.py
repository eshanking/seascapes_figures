from seascapes_figures.classes.population_class import Population
from seascapes_figures.utils import plotter, results_manager
import matplotlib.pyplot as plt
import numpy as np

p = Population(death_rate=0.0144,fitness_data='estimate')
fig,ax = plt.subplots(nrows=1,ncols=4,figsize = (8,4))

genotypes = np.array([0,1,2,3])
k = 0

for col in range(4):
    g = genotypes + k*4
    if k == 0:
        legend = True
    else:
        legend = False
    ax[col] = p.msw_grid(g,ax=ax[col],legend=legend)
    ax[col].set_xlabel(p.drug_units)
    k+=1

shift = 1
for col in range(1,4):
    ax[col] = p.shiftx(ax[col],0.01*shift)
    shift+=1

# ax[1,0].set_xlabel(p.drug_units)
# ax[1,1].set_xlabel(p.drug_units)

# ax[0].set_ylabel('Comparison',labelpad=35)

results_manager.save_fig(fig,'msw_grid_ecoli.pdf')

# fig = p.plot_msw(1,ncols=1,figsize=(10,2))
fig = p.plot_msw(1,ncols=1,figsize=(2.5,8),labelsize=14,ticklabelsize=12)
results_manager.save_fig(fig,'msw_1_ecoli.pdf')