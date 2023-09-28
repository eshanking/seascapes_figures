from fears.experiment import Experiment
import numpy as np

def make_data(plot=False,
              debug=False,
              n_sims=100):

    np.random.seed(2021)
    
    init_counts = np.zeros(16)
    init_counts[0] = 10**4


    options = {'doubling_time':.15,
               'death_rate':0.0144,
               'mut_rate':1.4*10**-8,
               'use_carrying_cap':True,
               'carrying_cap':10**9,
               'n_timestep':1700,
               'init_counts':init_counts,
               'k_abs':0.95,
               'k_elim':0.00839,
               'max_dose':50,
               'dose_schedule':24,
               'pad_right':True,
               'timestep_scale':1,
               'plot':plot,
               'dwell':True,
               'dwell_time':24,
               'regimen_length':21*24,
               'fitness_data':'from_file',
               'seascape_path':'results/seascape_library.xlsx',
               'plot_pop_size':True}
    
    p = np.array([0.0,0.2,0.4,0.6,0.8])

    n_sims = n_sims
    
    experiment_type = 'drug-regimen'
    
    e = Experiment(experiment_type=experiment_type,
                   n_sims=n_sims,
                   prob_drops=p,
                   population_options = options,
                   results_folder='results',
                   debug=debug)
    
    e.run_experiment()
    
    return e 

# if __name__ == '__main__':
#     make_data()