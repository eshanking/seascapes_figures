from scipy.stats.qmc import LatinHypercube
import scipy.interpolate
import numpy as np
import time
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from rate_survival_experiment_pharm import make_data
from fears.utils.stats import km_curve
# from seascapes_figures.utils import plotter

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

    km_data = km_curve(e)
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

def interpolate(x,y,grid=None,n=None,**kwargs):
    """Using scipy interpolate griddata, interpolates x,y over grid values.
    D dimensional data

    Args:
        x (length D tuple of 1-D ndarrays): data points
        y (ndarray): values associated with x data points
        grid (length D tuple of ndarrays): mesh grid to interpolate over
    """

    # if type(grid) is list:
    #     grid = tuple(grid)

    if grid is None:
        grid = []
        if n is None:
            n = 100
        for x_t in x:
            vect = np.linspace(min(x_t),max(x_t),num=n)
            grid.append(vect)
    
    grid = tuple(grid)
    grid = np.meshgrid(*grid)
    grid = tuple(grid)

    values = scipy.interpolate.griddata(x,y,grid,**kwargs)

    return values

sampler = LatinHypercube(d=num_dimensions)
sample = sampler.random(n=num_samples)
sample_raw = sample

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


###############################################################
# Plotting

# Interpolate results in 3D

x = tuple([sample_raw[:,0],sample_raw[:,1],sample_raw[:,2]])
v = interpolate(x,results,method='nearest')

# extent = np.min(x), np.max(x), np.min(y), np.max(y)

# death rate vs mutation rate

# sample = np.array([death_rate_sample,mut_rate_sample,carrying_cap_sample])

fig,ax = plt.subplots(nrows=1,ncols=3,constrained_layout=False,figsize=(10,15))
extent = np.min(x[0]), np.max(x[0]), np.min(x[1]), np.max(x[1])

ax[0].imshow(np.nanmean(v,axis=2),origin='lower')
# ax[0].title('Death rate vs mutation rate')
ax[0].set_xlabel('Death rate')
ax[0].set_ylabel('Mutation rate')

xt = ax[0].get_xticks()
xtl = scale_to_range(xt,death_rate_range)
ax[0].set_xticklabels(np.round(xtl,3))

yt = ax[0].get_yticks()
ytl = scale_to_range(yt,mut_rate_range)
ytl = [np.format_float_scientific(l,precision=1) for l in ytl]
ax[0].set_yticklabels(ytl)

# mutation rate vs carrying capacity

# extent = np.min(x[1]), np.max(x[1]), np.min(x[2]), np.max(x[2])

im = ax[1].imshow(np.nanmean(v,axis=0),origin='lower')
# ax[1].title('Mutation rate vs carrying capacity')
ax[1].set_xlabel('Mutation rate',labelpad=0.35)
ax[1].set_ylabel('Carrying capacity')

ax[1].set_xticklabels(ytl,rotation=45)

yt = ax[1].get_yticks()
ytl = scale_to_range(yt,carrying_cap_range)
ytl = [np.format_float_scientific(l,precision=1) for l in ytl]
ax[1].set_yticklabels(ytl)

carrying_cap_labels = ytl

# death rate vs carrying capacity

# extent = np.min(x[0]), np.max(x[0]), np.min(x[2]), np.max(x[2])

ax[2].imshow(np.nanmean(v,axis=1),origin='lower',interpolation='none')
# ax[1].title('Death rate vs carrying capacity')
ax[2].set_xlabel('Death rate')
ax[2].set_ylabel('Carrying capacity')

ax[2].set_yticklabels(carrying_cap_labels)

xt = ax[2].get_xticks()
xtl = scale_to_range(xt,death_rate_range)
ax[2].set_xticklabels(np.round(xtl,3))

# fig.colorbar(im,ax=ax[2],label='Range')
ax[1] = p.shiftx(ax[1],.1)
ax[2] = p.shiftx(ax[2],0.2)

axins = inset_axes(ax[1],
                    width="100%",  
                    height="5%",
                    loc='lower center',
                    borderpad=-7
                   )
fig.colorbar(im, cax=axins, orientation="horizontal",label='Survival probability range')

fig.savefig('lhs_analysis_nearest.pdf',bbox_inches='tight')