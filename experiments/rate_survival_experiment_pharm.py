from fears.experiment import Experiment
import numpy as np
    
def make_data(death_rate=None,
              mut_rate=None,
              n_sims=None,
              carrying_cap=None,
              debug=False,
              population_template=None,
              slopes=None):


   np.random.seed(2021)

   #  max_doses = [100]
   max_doses = [100]
   curve_types = ['pharm']
   experiment_type = 'rate-survival'
      
   if n_sims is None:
      n_sims = 1
   #  n_sims = 1

   # slopes = np.array([0.4,0.5,0.6,0.7])*10**-3

   # slopes = np.array([0.4])*10**-3

   #  slopes = np.array([0.5,1,2,4])*10**-4

   if death_rate is None:
      death_rate = 0.0144
   if mut_rate is None:
      mut_rate = 1.4*10**-8
   if carrying_cap is None:
      carrying_cap = 10**11
   if slopes is None:
      slopes = np.array([2,2.5,3,3.5])*10**-4

   init_counts = np.zeros(16)
   init_counts[0] = 10**10
    
   options = {'doubling_time':1,
            # 'death_rate':0.0144,
            'death_rate':death_rate,
         #    'mut_rate':10**-9,
            'mut_rate':mut_rate,
            # 'mut_rate':10**-3,
            'carrying_cap':True,
            'max_cells':carrying_cap,
            'n_timestep':2000,
            'init_counts':init_counts,
            # 'k_abs':0.95,
            # 'k_elim':0.00839,
            'k_elim':0,
         #    'max_dose':100,
            'dose_schedule':24,
            'pad_right':False,
            'timestep_scale':3,
            'plot':False,
         #    'ic50_data':'pyrimethamine_ic50.csv',
            'fitness_data':'estimate',
            'null_seascape':False,
            # 'null_seascape_dose':1
            # 'null_seascape_method':'sort'
            }


   e = Experiment(max_doses=max_doses,
                  slopes=slopes,
                  curve_types = curve_types,
                  experiment_type = experiment_type,
                  n_sims=n_sims,
                  passage = False,
                  passage_time = 96,
                  population_options=options,
                  population_template=population_template,
                  debug=debug)


   e.run_experiment()

   return e

if __name__ == '__main__':
    make_data()