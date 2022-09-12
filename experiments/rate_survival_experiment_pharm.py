from fears.experiment import Experiment
import numpy as np
    
def make_data(death_rate=None,
              mut_rate=None,
              n_sims=None,
              carrying_cap=None,
              debug=False,
              plot=False,
              population_template=None,
              timestep_scale=1,
              slopes=None):


   np.random.seed(2021)

   max_doses = [100]
   curve_types = ['pharm']
   experiment_type = 'rate-survival'
      
   if n_sims is None:
      n_sims = 1
   #  n_sims = 1

   if death_rate is None:
      death_rate = 0.0144
   if mut_rate is None:
      mut_rate = 1.4*10**-8
   if carrying_cap is None:
      carrying_cap = 10**11
   if slopes is None:
      slopes = np.array([12,16,20,24])*10**-3

   init_counts = np.zeros(16)
   init_counts[0] = 10**10
    
   options = {'doubling_time':.1,
            'death_rate':death_rate,
            'mut_rate':mut_rate,
            'use_carrying_cap':True,
            'carrying_cap':carrying_cap,
            'n_timestep':300,
            'init_counts':init_counts,
            'k_elim':0,
            # 'dose_schedule':24,
            'pad_right':False,
            'timestep_scale':timestep_scale,
            'plot':plot,
            'fitness_data':'from_file',
            'seascape_path':'results/seascape_library.xlsx',
            'null_seascape':False,
            'debug':False
            }

   e = Experiment(max_doses=max_doses,
                  slopes=slopes,
                  curve_types = curve_types,
                  experiment_type = experiment_type,
                  n_sims=n_sims,
                  passage = True,
                  passage_time = 24,
                  population_options=options,
                  population_template=population_template,
                  results_folder='results',
                  debug=debug)


   e.run_experiment()

   return e

# if __name__ == '__main__':
#     make_data()