from fears.population import Population
from fears.utils import plotter
import matplotlib.pyplot as plt
import numpy as np

np.random.seed(2022)

options = {'doubling_time':.15,
               'death_rate':0.0144,
               'mut_rate':1.4*10**-8,
               'use_carrying_cap':True,
               'carrying_cap':10**9,
               'n_timestep':1700,
            #    'init_counts':init_counts,
               'k_abs':0.95,
               'k_elim':0.00839,
               'max_dose':50,
               'dose_schedule':24,
               'pad_right':True,
               'timestep_scale':1,
            #    'plot':plot,
               'dwell':True,
               'dwell_time':24,
               'regimen_length':21*24,
               'fitness_data':'from_file',
               'seascape_path':'results/seascape_library.xlsx',
               'plot_pop_size':True}

p0 = Population(**options,curve_type='pulsed',prob_drop=0)

fig,ax_list = plt.subplots(nrows=2,figsize=(3,3))

ax = ax_list[0]
ax.plot(p0.drug_curve,color='black',linewidth=2)

ax = plotter.x_ticks_to_days(p0,ax)
ax.set_xlim(0,1000)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)

ax.set_ylabel('Concentration',fontsize=12)

ax_reg = ax_list[1]

x = np.arange(p0.n_timestep)

ax_reg.plot(x,p0.impulses,linewidth=0.5,color='black')
ax_reg = plotter.x_ticks_to_days(p0,ax_reg)
ax_reg.set_xlim(0,1000)

ax_reg = plotter.shrinky(ax_reg,0.3)
ax_reg = plotter.shifty(ax_reg,0.25)
ax_reg.spines["right"].set_visible(False)
ax_reg.spines["top"].set_visible(False)
ax_reg.spines["left"].set_visible(False)
ax_reg.set_yticks([])
ax_reg.set_xlabel('Days',fontsize=12)

fig.savefig('figures/schematic_drug_conc_curve.png',bbox_inches='tight')

drug_conc_range = [-4,4]
p1 = Population(fitness_data='random',
               n_allele=2,
               death_rate=0.1,
               drug_conc_range = drug_conc_range,
               ic50_limits=[-2.5,3],
               drugless_limits=[0.8,1.5])

p1.drugless_rates = [1.28949852, 1.14399848, 1.22802236, 0.93619847]
p1.ic50 = [-0.49205992, 1.76224515,  1.39341393,  2.84653598]

fig2,ax = plt.subplots(figsize=(2,2))
ax = plotter.plot_landscape(p,ax=ax,network_only=True)
fig2.savefig('figures/schematic_landscape.png',bbox_inches='tight')


fig3,ax = plt.subplots(figsize=(5,3))
fig3,ax = plotter.plot_fitness_curves(p,ax=ax,fig=fig3)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.legend(loc='best',frameon=False,fontsize=14)
fig3.savefig('figures/schematic_seascape.png',bbox_inches='tight')