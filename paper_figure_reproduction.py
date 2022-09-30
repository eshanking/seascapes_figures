from experiments import seascape_v_landscape, rate_survival_experiment_pharm, \
                        adherance_survival_experiment, multi_od_seascape_est, \
                        latin_hypercube_analysis
from figure_code import seascape_v_landscape_fig, ecoli_seascape, \
                        rate_of_change_km_fig, nonadherance_km_figure, \
                        example_timecourses, adh_histogram

"""Recreates all figures as they appear in the paper.

Note that this will take several hours to run and store a few gigabytes worth of data in 
the results folder.
"""

# Fig. 1

e_fig1 = seascape_v_landscape.make_data()
seascape_v_landscape_fig.make_fig(exp=e_fig1)

# Fig. 2

ecoli_seascape.make_fig()

# Fig. 3

roc_exp = rate_survival_experiment_pharm.make_data(n_sims=500)
adh_exp = adherance_survival_experiment.make_data(n_sims=500)
example_timecourses.make_fig(roc_exp=roc_exp,adh_exp=adh_exp)

# Fig. 4

rate_of_change_km_fig.make_fig(roc_exp=roc_exp)

# Fig. 5

nonadherance_km_figure.make_fig(adh_exp=adh_exp)

# Fig. 6

adh_histogram.make_fig(exp=adh_exp)

# Supplementals

multi_od_seascape_est.make_fig()
latin_hypercube_analysis.make_fig()