import sys
sys.path.append('/home/esk81')
from scipy.stats.qmc import LatinHypercube
import numpy as np
import time
import matplotlib.pyplot as plt
from seascapes_figures.experiments.rate_survival_experiment_pharm import make_data
from seascapes_figures.figure_code.rate_of_change_km_fig import make_fig

num_samples = 2
num_dimensions = 3 # death rate, mutation rate, and carrying capacity
n_sims = 10

per_sim_runtime = 14 # seconds
runtime_estimate = per_sim_runtime*num_samples*n_sims

print('Runtime estimate = ' + str(runtime_estimate))

def scale_to_range(x,range):
    width = max(range)-min(range)
    data_width = max(x) - min(x)

    x_scale = x/data_width
    x_scale = x_scale*width
    x_scale = x_scale + min(range)

    return x_scale

def get_outcome(e):

    result = make_fig(e)
    result = calculate_result_range(result,e)

    return result

def calculate_result_range(result,e):

    n_timestep = e.populations[0].n_timestep

    result = result['survival']
    # 
    num_conditions = len(result.keys())
    proportion_extinct = np.zeros(num_conditions)

    k = 0
    for key in result.keys():
        survival_times = result[key]
        num_extinct = len(np.argwhere(survival_times<n_timestep))
        proportion_extinct[k] = num_extinct/e.n_sims
    
    result_range = max(proportion_extinct) - min(proportion_extinct)

    return result_range

sampler = LatinHypercube(d=num_dimensions)
sample = sampler.random(n=num_samples)

# min_death_rate = 1/(12*24)
min_death_rate = 0.001
max_death_rate = 0.1
# max_death_rate = 1/(2*24)

death_rate_range = (min_death_rate,max_death_rate)

death_rate_sample = scale_to_range(sample[:,0],death_rate_range)

min_mut_rate = 10**-9
max_mut_rate = 10**-7

mut_rate_range = (min_mut_rate,max_mut_rate)
mut_rate_sample = scale_to_range(sample[:,1],mut_rate_range)

min_carrying_cap = 10**8
max_carrying_cap = 10**11

carrying_cap_range = (min_carrying_cap,max_carrying_cap)
carrying_cap_sample = scale_to_range(sample[:,2],carrying_cap_range)

carrying_cap_sample = [int(s) for s in carrying_cap_sample]

sample = np.array([death_rate_sample,mut_rate_sample,carrying_cap_sample])
sample = np.transpose(sample)

results = []

k = 0

tic = time.time()

for s in sample:
    death_rate = s[0]
    mut_rate = s[1]
    cc = s[2]

    e = make_data(death_rate=death_rate,mut_rate=mut_rate,n_sims=n_sims,carrying_cap=cc,debug=False)
    result = get_outcome(e)

    results.append(result)

    k+=1

toc = time.time()
elapsed = toc-tic
print(elapsed)

fig,ax = plt.subplots()

# mut rate v death rate
im = ax.scatter(sample[:,0],sample[:,1],c=results)

ax.scatter(0.014,1.4*10**-8,c='black',marker='x')

fig.colorbar(im,ax=ax,label='Range')

ax.set_ylabel('Mutation rate')
ax.set_xlabel('Death rate')

# ax.set_xticks([0.005,0.01,0.015,0.02])

fig.savefig('lhs_mutrate_vs_deathrate.pdf')

# death rate v carrying cap
fig,ax = plt.subplots()
im = ax.scatter(sample[:,0],sample[:,2],c=results)

ax.scatter(0.014,10**11,c='black',marker='x')

fig.colorbar(im,ax=ax,label='Range')

ax.set_ylabel('Carrying capacity')
ax.set_xlabel('Death rate')

# ax.set_xticks([0.005,0.01,0.015,0.02])

fig.savefig('lhs_carrying_cap_vs_death_rate.pdf')

# mut rate v carrying cap
fig,ax = plt.subplots()
im = ax.scatter(sample[:,1],sample[:,2],c=results)

ax.scatter(1.4*10**-8,10**11,c='black',marker='x')

fig.colorbar(im,ax=ax,label='Range')

ax.set_ylabel('Carrying capacity')
ax.set_xlabel('Mutation rate')

# ax.set_xticks([0.005,0.01,0.015,0.02])

fig.savefig('lhs_mutrate_vs_carrying_cap.pdf')