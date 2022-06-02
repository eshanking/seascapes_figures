from re import A
from seascapes_figures.classes.population_class import Population
from seascapes_figures.utils import results_manager
import matplotlib.pyplot as plt
import numpy as np

np.random.seed(2022)

p = Population(curve_type='pulsed',prob_drop=0.5,k_elim=0.02,k_abs=0.5,
                dose_schedule=12,max_dose=100)

fig,ax = plt.subplots()

ax.plot(p.drug_curve,color='black')
ax.set_xlim(0,500)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.set_ylabel('Drug concentration (uM)',fontsize=20)
ax.set_xlabel('Time (hr)',fontsize=20)

# xl = ax.get_xticklabels()
xt = ax.get_xticks()
ax.set_xticks(xt)
xt = [int(x) for x in xt]
ax.set_xticklabels(xt,fontsize=15)

yt = ax.get_yticks()
ax.set_yticks(yt)
yt = [int(x) for x in yt]
ax.set_yticklabels(yt,fontsize=15)

ax.set_ylim(0,600)

ax.annotate('$p_{forget} = 0.5$',(350,550),fontsize=18)

results_manager.save_fig(fig,'example_curve_p_forget_20.pdf')