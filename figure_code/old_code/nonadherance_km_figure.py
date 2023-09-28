import os
import matplotlib.pyplot as plt
import numpy as np
from fears.utils import results_manager, stats, plotter
import pandas as pd
import pickle


def make_fig(adh_exp=None,exp_info_path=None,resistance_outcome=[[1,2,4,8],[3,5,6,9,10,12]]):
    # suffix = '11182021_0001' # lab machine
    # suffix = '11112021_0000' # macbook
    
    # if suffix is None: 
    #     data_folder = adh_exp.results_path
    #     exp_info_file = adh_exp.experiment_info_path
    # else:
    #     data_folder = 'results_' + suffix
    #     exp_info_file = 'experiment_info_' + suffix + '.p'

    if adh_exp is None:
        adh_exp = pickle.load(open(exp_info_path,'rb'))
    
    fig,ax = plt.subplots(nrows=1,ncols=3,figsize=(8,2.5))
    
    # if suffix is None:
    #     exp_folders,exp_info = results_manager.get_experiment_results(exp=adh_exp)
    # else:
    #     exp_folders,exp_info = results_manager.get_experiment_results(suffix=suffix)
    exp_folders,exp_info = results_manager.get_experiment_results(exp=adh_exp)

    n_sims = exp_info.n_sims
    p_drop = exp_info.prob_drops
    
    exp_folders.reverse()
    p_drop = np.flip(p_drop)
    
    pop = exp_info.populations[0]
    tmax = int(pop.n_timestep)

    km_data = stats.km_curve(exp=adh_exp,resistance_outcome=resistance_outcome)

    #%%

    if type(resistance_outcome[0]) == list:
        key1 = 'resistance' + str(resistance_outcome[0])
    else:
        key1 = 'resistance ' + pop.int_to_binary(resistance_outcome[0])
    if type(resistance_outcome[1]) == list:
        key2 = 'resistance' + str(resistance_outcome[1])
    else:
        key2 = 'resistance ' + pop.int_to_binary(resistance_outcome[1])

    for k_abs in km_data.keys():
        
        exp_dict = km_data[k_abs]

        death_event_times = exp_dict['survival']
        gen1_resistance_times = exp_dict[key1]
        gen2_resistance_times = exp_dict[key2]


        ax[0] = plotter.plot_kaplan_meier(pop,death_event_times,
                                          ax=ax[0],
                                          n_sims=n_sims,
                                          label=k_abs,
                                          mode='survival',
                                          t_max=tmax)
        
        ax[1] = plotter.plot_kaplan_meier(pop,gen1_resistance_times,
                                          ax=ax[1],
                                          n_sims=n_sims,
                                          label=k_abs,
                                          mode='resistant',
                                          t_max=tmax)
        
        ax[2] = plotter.plot_kaplan_meier(pop,gen2_resistance_times,
                                          ax=ax[2],
                                          n_sims=n_sims,
                                        #   label=f"{k_abs_t:.1e}",
                                          label = k_abs,
                                          mode='resistant',
                                          t_max=tmax)
            
    
    for a in ax:
        a.spines["right"].set_visible(False)
        a.spines["top"].set_visible(False)
        a = plotter.x_ticks_to_days(pop,a)
        a.set_xlabel('Days')
        
    ax[2].legend(frameon=False,loc=[1.1,.3],title='$p_{forget}$',fontsize=8,ncol=1)
    
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
    ax[1].set_title('Single mutant',fontsize=8)
    ax[2].set_title('Double mutant',fontsize=8)

    # add vertical line to each axis

    x = np.ones(100)*21*24/pop.timestep_scale
    y = np.arange(100)    

    for a in ax:
        a.plot(x,y,'--',color='black',linewidth=2)

    # results_manager.save_fig(fig,'nonadherance_km_curve.pdf',bbox_inches='tight')
    fig.savefig('figures/adh_roc_curve.pdf',bbox_inches='tight')

    return km_data

#%%  perform pairwise log-rank tests and compute p values

    # analysis_keys = list(km_data.keys()) # endpoints being analyzed
    # experiment_keys = [str(p) for p in p_drop] # experiments performed
    
    # comparisons = [] # vector of all pairwise comparisons without duplicates
    
    # for i in range(len(experiment_keys)):
    #     j = i+1
    #     while j < len(experiment_keys):
    #         pair = (p_drop[i],p_drop[j])
    #         j+=1
    #         comparisons.append(pair)
    
    # p_values = {'survival':{},
    #             'resistance 0010':{},
    #             'resistance 0110':{}}
    
    # for ak in  analysis_keys:
    #     for pair in comparisons:
    #         key0 = str(pair[0])
    #         key1 = str(pair[1])
    #         sr = exp_info.log_rank_test(km_data[ak][key0],km_data[ak][key1])
    #         p_values[ak][str(pair)] = float(sr.p_value)#*n_tests # Mutliple hypothesis testing correction
    
    # p_values = pd.DataFrame(p_values)
    # result_path = dir_manager.make_resultspath_absolute(
    #     'nonadherance_km_curves_p_values.csv')
    
    # p_values.to_csv(result_path)
    
# if __name__ == '__main__':
#     make_fig() 

# def unpack(sim_path):

#     data_dict = pickle.load(open(sim_path,'rb'))
#     counts = data_dict['counts']
#     drug_conc = data_dict['drug_curve']
#     regimen = data_dict['regimen']

#     return counts, drug_conc, regimen