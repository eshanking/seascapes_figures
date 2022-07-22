import sys
sys.path.append('/home/esk81')
from scipy.stats.qmc import LatinHypercube
import numpy as np
import time
import matplotlib.pyplot as plt
from seascapes_figures.experiments.rate_survival_experiment_pharm import make_data
from seascapes_figures.figure_code.rate_of_change_km_fig import make_fig

np.random.seed(2022)
num_samples = 100
num_dimensions = 3 # death rate, mutation rate, and carrying capacity
n_sims = 10

slopes = [2*10**-4,3.5*10**-4]

per_sim_runtime = 4 # seconds
runtime_estimate = per_sim_runtime*num_samples*n_sims

print('Runtime estimate = ' + str(runtime_estimate))

def scale_to_range(data,range):

    scale_width = max(range)-min(range)
    data = np.array(data)
    # normalize data to between 0 and 1

    data_width = max(data) - min(data)

    data = data - min(data)

    data = data/data_width

    # scale normalized data to range
    data = data*scale_width
    data = data+min(range)

    return data

def get_outcome(e):

    km_data = make_fig(e)
    result = calculate_result_range(km_data,e)

    return result

def calculate_result_range(km_data,e):

    n_timestep = e.populations[0].n_timestep

    survival_km = km_data['survival']
    # 
    num_conditions = len(survival_km.keys())
    proportion_extinct = np.zeros(num_conditions)

    k = 0
    for key in survival_km.keys():
        survival_times = survival_km[key]
        num_extinct = len(np.argwhere(survival_times<n_timestep))
        proportion_extinct[k] = num_extinct/e.n_sims
        k+=1
    
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

# get a template population object by doing the experiment once
s =  sample[0]
death_rate = s[0]
mut_rate = s[1]
cc = s[2]

e = make_data(death_rate=death_rate,mut_rate=mut_rate,n_sims=n_sims,carrying_cap=cc,debug=False,slopes=slopes)

result = get_outcome(e)
results.append(result)

p = e.populations[0]

for s in sample[1:]:
    death_rate = s[0]
    mut_rate = s[1]
    cc = s[2]

    p.reset_drug_conc_curve(mut_rate=mut_rate,death_rate=death_rate,max_cells=cc)

    e = make_data(death_rate=death_rate,mut_rate=mut_rate,n_sims=n_sims,carrying_cap=cc,debug=False,
                  population_template=p,slopes=slopes)

    result = get_outcome(e)

    results.append(result)


toc = time.time()
elapsed = toc-tic
print(round(elapsed))

avg_runtime_per_sim = elapsed/(n_sims*num_samples)

print('Runtime per sim: ' + str(round(avg_runtime_per_sim)))

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