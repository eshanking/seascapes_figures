from seascapes_figures.classes.population_class import Population
from seascapes_figures.utils import plotter, results_manager
import matplotlib.pyplot as plt

p = Population()

fig,ax = plt.subplots(figsize=(5,5))

cbloc = [0.07,0.25,1,0.5]

ax = p.plot_landscape(ax=ax,
                      square=True,
                      node_size=800,
                      textsize=11,
                      cbloc=cbloc,
                      textcolor='white')

results_manager.save_fig(fig,'example_landscape.png')