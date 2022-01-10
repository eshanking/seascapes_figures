from seascapes_figures.classes.population_class import Population
from seascapes_figures.utils import plotter, results_manager
import matplotlib.pyplot as plt
import numpy as np

p = Population(death_rate=0.1)
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
        ax[row,col] = plotter.msw_grid(p,g,ax=ax[row,col],legend=legend)
        k+=1

plotter.shifty(ax[1,0],0.02)
plotter.shifty(ax[1,1],0.02)

ax[1,0].set_xlabel('Drug Concentration ($\u03BC$M)')
ax[1,1].set_xlabel('Drug Concentration ($\u03BC$M)')

ax[0,0].set_ylabel('Comparison',labelpad=20)
ax[1,0].set_ylabel('Comparison',labelpad=20)

results_manager.save_fig(fig,'msw_grid.pdf')

fig = plotter.plot_msw(p,4,ncols=1)
results_manager.save_fig(fig,'msw_4.pdf')