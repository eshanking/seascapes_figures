from seascapes_figures.classes.population_class import Population
from seascapes_figures.utils import plotter, results_manager
import matplotlib.pyplot as plt

p = Population()

fig,ax = plt.subplots()

ax = p.plot_landscape(p,network_only=True,
                            ax=ax,
                            square=True,
                            node_size=400,
                            textsize=7,
                            plot_sub_network=True,
                            sub_network=[1,3,5,9],
                            sub_network_color='black')

results_manager.save_fig(fig,'neighbors_genotype_1.pdf')