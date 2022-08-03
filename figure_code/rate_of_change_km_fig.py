import os
import matplotlib.pyplot as plt
import numpy as np
from fears.utils import results_manager, dir_manager, stats
import pandas as pd
import pickle

def make_fig(roc_exp=None,exp_info_path=None,save=True,resistance_outcome=[14,15]):

    if roc_exp is None:
        roc_exp = pickle.load(open(exp_info_path,'rb'))

    data_folder = roc_exp.results_path
    exp_info_file = roc_exp.experiment_info_path

    exp_folders,exp_info = results_manager.get_experiment_results(exp=roc_exp)

    n_sims = exp_info.n_sims
    k_abs = exp_info.slopes
    
    fig,ax = plt.subplots(nrows=1,ncols=3,figsize=(8,2.5))
    
    pop = exp_info.populations[0]

    km_data = stats.km_curve(exp=roc_exp,resistance_outcome=resistance_outcome)
    tmax = int(pop.n_timestep)

    key1 = 'resistance ' + pop.int_to_binary(resistance_outcome[0])
    key2 = 'resistance ' + pop.int_to_binary(resistance_outcome[1])

    for k_abs in km_data.keys():
        
        exp_dict = km_data[k_abs]

        death_event_times = exp_dict['survival']
        gen1_resistance_times = exp_dict[key1]
        gen2_resistance_times = exp_dict[key2]


        ax[0] = pop.plot_kaplan_meier(death_event_times,
                                          ax=ax[0],
                                          n_sims=n_sims,
                                          label=k_abs,
                                          mode='survival',
                                          t_max=tmax)
        
        ax[1] = pop.plot_kaplan_meier(gen1_resistance_times,
                                          ax=ax[1],
                                          n_sims=n_sims,
                                          label=k_abs,
                                          mode='resistant',
                                          t_max=tmax)
        
        ax[2] = pop.plot_kaplan_meier(gen2_resistance_times,
                                          ax=ax[2],
                                          n_sims=n_sims,
                                        #   label=f"{k_abs_t:.1e}",
                                          label = k_abs,
                                          mode='resistant',
                                          t_max=tmax)
    
    for a in ax:
        a.spines["right"].set_visible(False)
        a.spines["top"].set_visible(False)
        
    # ax[2].legend(frameon=False,loc=[-2,-0.5],title='$k_{abs}$',fontsize=8,ncol=4)
    ax[2].legend(frameon=False,loc=[1.1,.3],title='$k_{abs}$',fontsize=8,ncol=1)
    
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
    title_t = 'Resistant genotype = ' + pop.int_to_binary(resistance_outcome[0])
    ax[1].set_title(title_t,fontsize=8)

    title_t = 'Resistant genotype = ' + pop.int_to_binary(resistance_outcome[1])
    ax[2].set_title(title_t,fontsize=8)
    
    # max sure all x lims are the same
    xmax = ax[0].get_xlim()[1]
    for a in ax:
        if a.get_xlim()[1] > xmax:
            xmax = a.get_xlim()[1]
    
    xmax = xmax/2

    for a in ax:
        a.set_ylim([-10,110])
        # a.set_xlim([0,xmax])
        a = pop.x_ticks_to_days(a)

    if save:
        results_manager.save_fig(fig,'roc_km_curve.pdf',bbox_inches='tight')

    return fig,ax
    #%%
    # perform pairwise log-rank tests and compute p values

    # if save:
    #     analysis_keys = list(km_data.keys()) # endpoints being analyzed
    #     experiment_keys = [str(p) for p in k_abs] # experiments performed
        
    #     comparisons = [] # vector of all pairwise comparisons without duplicates
        
    #     for i in range(len(experiment_keys)):
    #         j = i+1
    #         while j < len(experiment_keys):
    #             pair = (k_abs[i],k_abs[j])
    #             j+=1
    #             comparisons.append(pair)
        
    #     p_values = {'survival':{},
    #                 'resistance 0111':{},
    #                 'resistance 1111':{}}
        
    #     n_tests = len(k_abs)-1
        
    #     for ak in  analysis_keys:
    #         for pair in comparisons:
    #             key0 = str(pair[0])
    #             key1 = str(pair[1])
    #             sr = exp_info.log_rank_test(km_data[ak][key0],km_data[ak][key1])
    #             p_values[ak][str(pair)] = float(sr.p_value)*n_tests # Mutliple hypothesis testing correction
        
    #     p_values = pd.DataFrame(p_values)
    #     result_path = dir_manager.make_resultspath_absolute(
    #         'rate_of_change_km_curves_p_values.csv')
    
    
    # p_values.to_csv(result_path)

    # return km_data
# #%%
# if __name__ == '__main__':

#     eip = '/Users/kinge2/Library/CloudStorage/Box-Box/seascapes_figures/results/results_06222022_0000/experiment_info_06222022_0000.p'
#     # eip = '/Users/kinge2/Library/CloudStorage/Box-Box/seascapes_figures/results/results_06222022_0001/experiment_info_06222022_0001.p'
#     make_fig(exp_info_path=eip)
# #%%     
