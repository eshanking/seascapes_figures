from scipy.stats.qmc import LatinHypercube
import numpy as np
from seascapes_figures.experiments.rate_survival_experiment_pharm import make_data
from seascapes_figures.figure_code.rate_of_change_km_fig import make_fig

num_samples = 2

def scale_to_range(x,range):
    width = max(range)-min(range)
    data_width = max(x) - min(x)

    x_scale = x/data_width
    x_scale = x_scale*width
    x_scale = x_scale + min(range)

    return x_scale

def get_outcome(e):

    result = make_fig(e)

    return result

sampler = LatinHypercube(d=2)
sample = sampler.random(n=num_samples)

min_death_rate = 1/(12*24)
max_death_rate = 1/(2*24)

death_rate_range = (min_death_rate,max_death_rate)

death_rate_sample = scale_to_range(sample[:,0],death_rate_range)

min_mut_rate = 10**-9
max_mut_rate = 10**-7

mut_rate_range = (min_mut_rate,max_mut_rate)
mut_rate_sample = scale_to_range(sample[:,1],mut_rate_range)

sample = np.array([death_rate_sample,mut_rate_sample])
sample = np.transpose(sample)

results = np.zeros(num_samples)

k = 0
for s in sample:
    death_rate = s[0]
    mut_rate = s[1]

    e = make_data(death_rate=death_rate,mut_rate=mut_rate,n_sims=100)
    result = get_outcome(e)
    results[k] = result
    k+=1
    