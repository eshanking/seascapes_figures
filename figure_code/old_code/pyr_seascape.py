from seascapes_figures.utils import plotter, results_manager
from seascapes_figures.classes.population_class import Population
import matplotlib.pyplot as plt

def make_figure():
    fig,ax = plt.subplots(nrows=1,ncols=1,sharey=False,figsize=(7,4))
    
    p_pyr = Population(ic50_data='pyrimethamine_ic50.csv')
    # p_cyc = Population(ic50_data='cycloguanil_ic50.csv')
    
    t,ax = plotter.plot_fitness_curves(p_pyr,ax=ax,show_legend=False,
                                        linewidth=2)
    # t,ax[1] = plotter.plot_fitness_curves(p_cyc,ax=ax[1],show_legend=False,
    #                                     linewidth=2)

    conc = [10**-3, 10**1, 10**5]

    for c in conc:

        ax,lax = plotter.add_landscape_to_fitness_curve(c,ax,p_pyr,
                                                        position='bottom',
                                                        width=4.5,
                                                        pad=2.2,
                                                        node_size=500,
                                                        textsize=8,
                                                        square=True,
                                                        colorbar=False,
                                                        textcolor='white')
    
    # ax,lax = plotter.add_landscape_to_fitness_curve(conc[2],ax,p_pyr,
    #                                             position='bottom',
    #                                             width=4.5,
    #                                             pad=2.1,
    #                                             node_size=500,
    #                                             textsize=9,
    #                                             square=True,
    #                                             colorbar=False)

    xt = [-3,-1,1,3,5]
    xt = [10**p for p in xt]
    ax.set_xticks(xt)
    ax.xaxis.tick_top()
    ax.xaxis.set_label_coords(0.5, 1.25)
    ax.set_ylabel('Growth rate ($hr^{-1}$)')

    
    ax.legend(ncol=2,frameon=False,loc=(1,0.3),title='Genotypes')
    
    results_manager.save_fig(fig,'pyr_seascapes.pdf')

if __name__ == '__main__':
    make_figure()