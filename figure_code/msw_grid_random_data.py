from fears.population import Population
from fears.utils import plotter, dir_manager
import numpy as np
import matplotlib.pyplot as plt

np.random.seed(109)
drug_conc_range = [-3,5]
p = Population(fitness_data='random',
               n_allele=2,
               death_rate=0.1,
               drug_conc_range = drug_conc_range,
               ic50_limits=[-6.5,-3],
               drugless_limits=[0.8,1.5])

fig,ax = plt.subplots()

genotypes = np.array([0,1,2,3])
ax = plotter.msw_grid(p,genotypes,ax=ax,
                      labelsize=12,
                      ticklabelsize=12,
                      comp_annotate_pos=10**-3.5,
                      legendloc=(0,1.05))

fig.savefig('2_allele_grid.pdf')

# fig2 = plotter.plot_msw(p,2,ncols=1,figsize=(2.5,4))
fig2 = plotter.plot_msw(p,2,ncols=2,figsize=(6,2))
fig2.savefig('2_allele_msws.pdf')

fig3,ax = plt.subplots(figsize=(6,3))

fig3,ax = plotter.plot_fitness_curves(p,ax=ax,fig=fig3)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.legend(loc='best',frameon=False,fontsize=14)
fig3.savefig('random_seascape.pdf')