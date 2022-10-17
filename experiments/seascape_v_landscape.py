from fears.experiment import Experiment
import numpy as np

def make_data():
    np.random.seed(109)
    
    init_counts = np.zeros(4)
    init_counts[0] = 10**10
    
    transition_times = [1000,2000]

    ic50 = [3.64972198, 3.24736065, 1.19048008, 4.8325847 ]
    drugless_rates = [1.30573026, 1.2457132 , 1.3496418 , 1.09728462]

    options = {'doubling_time':1.5,
               'death_rate':0,
               'mut_rate':10**-9,
               'n_sims':10,
               # 'carrying_cap':False,
               'max_cells':10**11,
               'n_timestep':3000,
               'init_counts':init_counts,
               'timestep_scale':2,
               'plot':False,
               'plot_drug_curve':False,
               'drug_log_scale':True,
               'counts_log_scale':True,
               'max_dose':10**0,
               'drugless_rates':drugless_rates,
               'ic50':ic50,
               'drug_conc_range':(-3,5)}
    
    e = Experiment(debug=False,
                   experiment_type='ramp_up_down',
                   n_allele=2,
                   ramp = 1,
                   transition_times=transition_times,
                   population_options=options,
                   # null_seascape_dose=10**-3,
                   null_seascape_dose=10**1,
                   first_dose=10**-3,
                   second_dose=10**1,
                   third_dose=10**5)
    
    e.run_experiment()
    
    return e

# if __name__ == '__main__':
#     make_data()