import os
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoLocator
import numpy as np
from seascapes_figures.utils import results_manager, plotter, dir_manager
import pandas as pd
import pickle

def make_fig(roc_exp=None,exp_info_path=None):

    if roc_exp is None:
        roc_exp = pickle.load(open(exp_info_path,'rb'))

    data_folder = roc_exp.results_path
    exp_info_file = roc_exp.experiment_info_path
    
    # exp_folders,exp_info = results_manager.get_experiment_results(data_folder,
    #                                                              exp_info_file)

    exp_folders,exp_info = results_manager.get_experiment_results(exp=roc_exp)

    max_cells = exp_info.populations[0].max_cells
    n_sims = exp_info.n_sims
    k_abs = exp_info.slopes
    
    fig,ax = plt.subplots(nrows=1,ncols=3,figsize=(8,2.5))
    
    pop = exp_info.populations[0]
    
    km_data = {'survival':{},
               'resistance 0111':{},
               'resistance 1111':{}}
    
    for exp in exp_folders:
    
        k_abs_t = exp[exp.find('=')+1:]
        k_abs_t = k_abs_t.replace(',','.')
        k_abs_t = float(k_abs_t)
        
        # print(f"{k_abs_t:.2e}")

        num = np.argwhere(k_abs == k_abs_t)
        num = num[0,0]
        
        sim_files = os.listdir(path=exp)
        sim_files = sorted(sim_files)
        
        # KM data 
        death_event_obs = np.zeros(n_sims)
        death_event_times = np.zeros(n_sims)
        
        # time to genotype 2
        # gen14_resistance_obs = np.zeros(n_sims)
        # gen14_resistance_times = np.zeros(n_sims)
        gen7_resistance_obs = np.zeros(n_sims)
        gen7_resistance_times = np.zeros(n_sims)
        
        # time to genotype 6
        gen15_resistance_obs = np.zeros(n_sims)
        gen15_resistance_times = np.zeros(n_sims)
        
        k=0
        while k < len(sim_files):
        # while k < 10:
            sim = sim_files[k]
            sim = exp + os.sep + sim
            # print(sim)
            data_dict = results_manager.get_data(sim)
            # dc = data[:,-1]
            # data = data[:,0:-1]

            dc = data_dict['drug_curve']
            data = data_dict['counts']
            
            death_event_obs[k],death_event_times[k] = \
                exp_info.extinction_time(pop,data,thresh=1)
                
            # gen14_resistance_obs[k],gen14_resistance_times[k] = \
            #     exp_info.resistance_time(pop,data,14,thresh=.1)

            gen7_resistance_obs[k],gen7_resistance_times[k] = \
                exp_info.resistance_time(pop,data,7,thresh=.1)
    
            gen15_resistance_obs[k],gen15_resistance_times[k] = \
                exp_info.resistance_time(pop,data,15,thresh=.1)
                
            k+=1

        # get lenth of analysis time
        tmax = int(pop.n_timestep)

        ax[0] = pop.plot_kaplan_meier(death_event_times,
                                          ax=ax[0],
                                          n_sims=n_sims,
                                          label=str(k_abs_t),
                                          mode='survival',
                                          t_max=tmax)
        
        ax[1] = pop.plot_kaplan_meier(gen7_resistance_times,
                                          ax=ax[1],
                                          n_sims=n_sims,
                                          label=str(k_abs_t),
                                          mode='resistant',
                                          t_max=tmax)
        
        ax[2] = pop.plot_kaplan_meier(gen15_resistance_times,
                                          ax=ax[2],
                                          n_sims=n_sims,
                                          label=f"{k_abs_t:.1e}",
                                          mode='resistant',
                                          t_max=tmax)
        
        km_data['survival'][str(k_abs_t)] = death_event_times
        km_data['resistance 0111'][str(k_abs_t)] = gen7_resistance_times
        km_data['resistance 1111'][str(k_abs_t)] = gen15_resistance_times
    
    for a in ax:
        a.spines["right"].set_visible(False)
        a.spines["top"].set_visible(False)
        
    ax[2].legend(frameon=False,loc='upper left',title='$k_{abs}$',fontsize=8)
    
    pad = 0.05
    
    ax[0].set_ylabel('% surviving')
    pos1 = ax[1].get_position()
    pos1.x0 = pos1.x0 + pad
    pos1.x1 = pos1.x1 + pad
    ax[1].set_position(pos1)
    
    pos2 = ax[2].get_position()
    pos2.x0 = pos2.x0 + pad*2
    pos2.x1 = pos2.x1 + pad*2
    ax[2].set_position(pos2)
    
    ax[0].set_title('Survival of infectious agent',fontsize=8)
    ax[1].set_title('Resistant genotype = 0111',fontsize=8)
    ax[2].set_title('Resistant genotype = 1111',fontsize=8)
    
    # max sure all x lims are the same
    xmax = ax[0].get_xlim()[1]
    for a in ax:
        if a.get_xlim()[1] > xmax:
            xmax = a.get_xlim()[1]
    
    xmax = xmax/2

    for a in ax:
        a.set_ylim([-10,110])
        a.set_xlim([0,xmax])
        a = pop.x_ticks_to_days(a)

    results_manager.save_fig(fig,'roc_km_curve.pdf',bbox_inches='tight')
    #%%
    # perform pairwise log-rank tests and compute p values
    analysis_keys = list(km_data.keys()) # endpoints being analyzed
    experiment_keys = [str(p) for p in k_abs] # experiments performed
    
    comparisons = [] # vector of all pairwise comparisons without duplicates
    
    for i in range(len(experiment_keys)):
        j = i+1
        while j < len(experiment_keys):
            pair = (k_abs[i],k_abs[j])
            j+=1
            comparisons.append(pair)
    
    p_values = {'survival':{},
                'resistance 0111':{},
                'resistance 1111':{}}
    
    n_tests = len(k_abs)-1
    
    for ak in  analysis_keys:
        for pair in comparisons:
            key0 = str(pair[0])
            key1 = str(pair[1])
            sr = exp_info.log_rank_test(km_data[ak][key0],km_data[ak][key1])
            p_values[ak][str(pair)] = float(sr.p_value)*n_tests # Mutliple hypothesis testing correction
    
    p_values = pd.DataFrame(p_values)
    result_path = dir_manager.make_resultspath_absolute(
        'rate_of_change_km_curves_p_values.csv')
    
    p_values.to_csv(result_path)
#%%
if __name__ == '__main__':

    eip = '/Users/kinge2/repos/seascapes_figures/results/results_01212022_0003/experiment_info_01212022_0003.p'    
    make_fig(exp_info_path=eip)
#%%     