import sys
sys.path.append('/Users/eshanking/repos')
from seascapes_figures.classes.population_class import Population
from seascapes_figures.utils import plotter, results_manager
import matplotlib.pyplot as plt
import numpy as np

p = Population(death_rate=0.01,fitness_data='estimate')
fig,ax = plt.subplots(nrows=2,ncols=2,figsize = (6,8))

genotypes = np.array([0,1,2,3])
k = 0

for row in range(2):
    for col in range(2):
        g = genotypes + k*4
        if k == 0:
            legend = True
        else:
            legend = False
        ax[row,col] = p.msw_grid(g,ax=ax[row,col],legend=legend)
        k+=1

p.shifty(ax[1,0],0.02)
p.shifty(ax[1,1],0.02)

ax[1,0].set_xlabel(p.drug_units)
ax[1,1].set_xlabel(p.drug_units)

ax[0,0].set_ylabel('Comparison',labelpad=20)
ax[1,0].set_ylabel('Comparison',labelpad=20)

results_manager.save_fig(fig,'msw_grid_ecoli.pdf')

fig = p.plot_msw(4,ncols=1)
results_manager.save_fig(fig,'msw_4_ecoli.pdf')