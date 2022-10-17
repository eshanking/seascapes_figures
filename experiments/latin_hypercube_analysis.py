#%%
from scipy.stats.qmc import LatinHypercube, scale
from scipy import stats
import numpy as np
import pandas as pd
import time
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from fears.experiment import Experiment
from fears.utils import plotter

#%%

def make_fig():
    np.random.seed(2022)
    num_samples = 1000
    num_dimensions = 3 # death rate, mutation rate, and carrying capacity
    n_sims = 10

    # slopes = [2*10**-4,3.5*10**-4]
    slopes = [0.001]

    per_sim_runtime = 0.5 # seconds
    runtime_estimate = per_sim_runtime*num_samples*n_sims

    print('Runtime estimate = ' + str(runtime_estimate))

    sampler = LatinHypercube(d=num_dimensions)
    sample = sampler.random(n=num_samples)
    sample_raw = sample

    # min_death_rate = 1/(12*24)
    min_death_rate = 0.001
    max_death_rate = 0.1
    # # max_death_rate = 1/(2*24)

    min_mut_rate = -9
    max_mut_rate = -7

    min_carrying_cap = 8
    max_carrying_cap = 11

    lbounds = [min_death_rate,min_mut_rate,min_carrying_cap]
    ubounds = [max_death_rate,max_mut_rate,max_carrying_cap]

    sample = scale(sample,lbounds,ubounds)

    sample[:,1] = [10**s for s in sample[:,1]] # mutation rate from log to linear  
    sample[:,2] = [int(10**s) for s in sample[:,2]] # carrying capacity from log to linear

    results = []

    tic = time.time()

    # get a template population object by doing the experiment once
    s =  sample[0]
    death_rate = s[0]
    mut_rate = s[1]
    cc = s[2]

    init_counts = np.zeros(16)
    init_counts[0] = cc/100

    options = {'doubling_time':.15,
            'death_rate':death_rate,
            'mut_rate':mut_rate,
            'use_carrying_cap':True,
            'carrying_cap':cc,
            'n_timestep':500,
            'init_counts':init_counts,
            'k_elim':0,
            'fitness_data':'from_file',
            'seascape_path':'results/seascape_library.xlsx',
            'plot':False
            }

    e = Experiment(max_doses=[100],
                    slopes=slopes,
                    curve_types = ['pharm'],
                    experiment_type = 'rate_survival_lhs',
                    n_sims=n_sims,
                    passage = True,
                    passage_time = 24,
                    population_options=options,
                    results_folder='results',
                    debug=False)

    res = e.run_experiment()
    results.append(res)

    #%%
    p = e.populations[0]

    cnt = 0

    for s in sample[1:]:

        death_rate = s[0]
        mut_rate = s[1]
        cc = s[2]

        init_counts[0] = cc/100

        p.reset_drug_conc_curve(mut_rate=mut_rate,death_rate=death_rate,
                                max_cells=cc,init_counts=init_counts)

        e = Experiment(max_doses=[100],
                    slopes=slopes,
                    curve_types = ['pharm'],
                    experiment_type = 'rate_survival_lhs',
                    n_sims=n_sims,
                    passage = True,
                    passage_time = 24,
                    population_options=options,
                    results_folder='results',
                    population_template=p,
                    debug=False)

        res = e.run_experiment()

        results.append(res)
        
        pct = int(100*cnt/num_samples)
        if np.mod(cnt,5) == 0:
            print('Progress = ' + str(pct) + '%')
        cnt+=1

    toc = time.time()
    elapsed = toc-tic
    # print(round(elapsed))

    avg_runtime_per_sim = elapsed/(n_sims*num_samples)

    print('Runtime per sim: ' + str(round(avg_runtime_per_sim)))

    #%% Save data

    data_dict = {}
    data_dict['death rate'] = sample[:,0]
    data_dict['mutation rate'] = sample[:,1]
    data_dict['carrying capacity'] = sample[:,2]
    data_dict['result'] = results

    df = pd.DataFrame(data_dict)

    date_str = time.strftime('%m%d%Y',time.localtime())

    df.to_csv('results/lhs_analysis_' + date_str + '.csv')

    #%% Statistics

    # corr = df.corr()

    # corr = corr['result']

    # corr_array = corr.values

    # t_array = [r*(((num_samples-2)/(1-r**2))**0.5) for r in corr_array[:-1]]

    # dof = num_samples-2 # degrees of freedom

    # p = [stats.t.sf(abs(t),dof) for t in t_array]

    # p.append(1)

    # stats_df = pd.DataFrame(p,index=corr.index)

    # stats_df = pd.concat([stats_df,corr],axis=1)

    # stats_df.columns = ['p values','correlation']

    # stats_df.to_csv('results/lhs_analysis_correlations_' + date_str + '.csv')
    #%% Plotting

    fig,ax = plt.subplots(nrows=1,ncols=3,constrained_layout=False,figsize=(8,3))

    # High death rate
    results = np.array(results)
    vmin = min(results)
    vmax = max(results)

    death_rate_low = (max_death_rate-min_death_rate)/3 + min_death_rate

    # Filter by death rate
    death_rate_sample = sample[:,0]

    low_indx = death_rate_sample <= death_rate_low

    # s_low = sample[low_indx,:]
    res_low = results[low_indx]

    ax[0].scatter(sample[low_indx,1],sample[low_indx,2],c=results[low_indx],vmin=vmin,vmax = vmax)
    ax[0].set_title('Low death rate',fontsize=12)


    death_rate_med = (max_death_rate-min_death_rate)/3 + death_rate_low

    med_indx = death_rate_sample <= death_rate_med
    med_indx[low_indx] = False

    ax[1].scatter(sample[med_indx,1],sample[med_indx,2],c=results[med_indx],vmin=vmin,vmax = vmax)
    ax[1].set_title('Medium death rate',fontsize=12)

    high_indx = death_rate_sample > death_rate_med


    sc = ax[2].scatter(sample[high_indx,1],sample[high_indx,2],c=results[high_indx],vmin=vmin,vmax = vmax)
    ax[2].set_title('High death rate',fontsize=12)

    for a in ax:
            # a = ax[r,c]
        a.set_xlim(10**min_mut_rate,10**max_mut_rate)
        a.set_ylim(10**min_carrying_cap,10**max_carrying_cap)
        a.tick_params(axis='both', labelsize=12)
        a.set_xscale('log')
        a.set_yscale('log')
        a.set_ylabel('Carrying capacity',fontsize=12)
        a.set_xlabel('Mutation rate',fontsize=12)

    axins = inset_axes(ax[2],
                        width="10%",  
                        height="100%",
                        loc='right',
                        borderpad=-3
                    )

    fig.colorbar(sc, cax=axins, orientation="vertical",label='Survival probability')

    fig.tight_layout()

    fig.savefig('figures/lhs_analysis_' + date_str + '.pdf',bbox_inches='tight')