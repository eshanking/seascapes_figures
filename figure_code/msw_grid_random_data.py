from fears.population import Population
from fears.utils import plotter
import numpy as np
import matplotlib.pyplot as plt

np.random.seed(109)
drug_conc_range = [-3,5]
p = Population(fitness_data='random',
               n_allele=2,
               death_rate=0.5,
               drug_conc_range = drug_conc_range,
               ic50_limits=[-6.5,-3],
               drugless_limits=[0.8,1.5])

p.plot_fitness_curves()

fig,ax = plt.subplots()

genotypes = np.array([0,1,2,3])
ax = plotter.msw_grid(p,genotypes,ax=ax,
                      labelsize=12,
                      ticklabelsize=12,
                      comp_annotate_pos=10**-3.5,
                      legendloc=(0,1.05))

