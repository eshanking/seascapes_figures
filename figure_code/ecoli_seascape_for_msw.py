from seascapes_figures.classes.population_class import Population
from seascapes_figures.utils import results_manager
import matplotlib.pyplot as plt
import numpy as np

p = Population(fitness_data='estimate')

fig,ax = plt.subplots(figsize=(6,3))

fig,ax = p.plot_fitness_curves()

vert_lines_ydata = np.arange(9)/10
vert_lines_kwargs = {'linewidth':3,'alpha':0.7}

width = 2.5

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.set_ylabel('Growth rate ($hr^{-1}$)',fontsize=20)
ax.set_xlabel('Drug concentration (μg/mL)',fontsize=20)

ax.set_xlim([10**-2.7,10**3])
results_manager.save_fig(fig,'ecoli_seascape_for_msw.pdf')