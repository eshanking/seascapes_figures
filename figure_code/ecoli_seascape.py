from seascapes_figures.classes.population_class import Population
from seascapes_figures.utils import results_manager
import matplotlib.pyplot as plt
import numpy as np

p = Population(fitness_data='estimate')

fig,ax = plt.subplots(figsize=(6,4))

fig,ax = p.plot_fitness_curves()

vert_lines_ydata = np.arange(9)/10
vert_lines_kwargs = {'linewidth':3,'alpha':0.7}

ax,lax = p.add_landscape_to_fitness_curve(1,ax,
                                          width=3.5,
                                          height=0.7,
                                          pad=-.25,
                                          colorbar=False,
                                          vert_lines_ydata=vert_lines_ydata,
                                          vert_lines_kwargs=vert_lines_kwargs)

ax,lax = p.add_landscape_to_fitness_curve(10**-3,ax,
                                          width=3.5,
                                          height=0.7,
                                          pad=-.25,
                                          colorbar=False,
                                          vert_lines_ydata=vert_lines_ydata,
                                          vert_lines_kwargs=vert_lines_kwargs)

ax,lax = p.add_landscape_to_fitness_curve(10**3,ax,
                                          width=3.5,
                                          height=0.7,
                                          pad=-.25,
                                          colorbar=True,
                                          cbloc = [0.55,0.25,0.5,0.4],
                                          vert_lines_ydata=vert_lines_ydata,
                                          vert_lines_kwargs=vert_lines_kwargs)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.set_ylabel('Growth rate ($hr^{-1}$)',fontsize=20)
ax.set_xlabel('Drug concentration (μg/mL)',fontsize=20)

results_manager.save_fig(fig,'ecoli_seascape_with_landscapes.pdf')