from fears.population import Population
from fears.utils import plotter
import matplotlib.pyplot as plt
import numpy as np

def make_fig():
    p = Population(fitness_data='from_file',seascape_path='results/seascape_library.xlsx')

    fig,ax = plt.subplots(figsize=(11,5.5))

    xdata = np.logspace(-3.5,4.5,num=10000)
    fig,ax = p.plot_fitness_curves(fig=fig,ax=ax,legend_cols=2,xdata=xdata);

    vert_lines_ydata = np.arange(9)/10
    vert_lines_kwargs = {'linewidth':3,'alpha':0.7}

    width = 2.5

    ax,lax = plotter.add_landscape_to_fitness_curve(10**0,ax,p,
                                            width=width,
                                            height=0.1,
                                            pad=0,
                                            ypos=0.11,
                                            colorbar=False,
                                            vert_lines_ydata=vert_lines_ydata,
                                            vert_lines_kwargs=vert_lines_kwargs,
                                            textcolor='w')

    ax,lax = plotter.add_landscape_to_fitness_curve(10**-2.5,ax,p,
                                            width=width,
                                            height=0.1,
                                            pad=0,
                                            ypos=0.11,
                                            colorbar=False,
                                            vert_lines_ydata=vert_lines_ydata,
                                            vert_lines_kwargs=vert_lines_kwargs,
                                            textcolor='w')

    ax,lax = plotter.add_landscape_to_fitness_curve(10**2.5,ax,p,
                                            width=width,
                                            height=0.1,
                                            pad=0,
                                            ypos=0.11,
                                            colorbar=True,
                                            cbloc = [0.7,0.25,1,0.5],
                                            vert_lines_ydata=vert_lines_ydata,
                                            vert_lines_kwargs=vert_lines_kwargs,
                                            textcolor='w')

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.set_ylabel('Growth rate ($hr^{-1}$)',fontsize=20)
    ax.set_xlabel('Drug concentration (μg/mL)',fontsize=20)

    ax.set_ylim(0,0.12)
    ax.set_xlim([10**-2.7,10**3.5])

    fig.savefig('figures/seascape_with_landscapes.pdf',bbox_inches='tight')
    return p
# results_manager.save_fig(fig,'ecoli_seascape_with_landscapes.pdf')
make_fig()