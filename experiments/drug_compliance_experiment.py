from seascapes_figures.classes.experiment_class import Experiment
import numpy as np
# import matplotlib.pyplot as plt

def make_data():
    np.random.seed(2022)

    options = {'n_impulse':20,
            'k_abs':0.04,
            'k_elim':0.03,
            'max_dose':150,
            'n_timestep':1000,
            'mut_rate':1.4*10**-8,
            'death_rate':0.0144}

    e = Experiment(experiment_type='drug-regimen',
                    n_sims=10,
                    prob_drops=[0,.4,.7],
                    population_options=options,
                    debug=True)

    e.run_experiment()
    e.plot_barchart()
    e.save_images()

    # import matplotlib.pyplot as plt

    # fig,ax = plt.subplots()
    # ax.bar([0,1,2],[10,20,30])
    # ax.set_xticks([0,1,2])
    # ax.set_xticklabels(['a','b','c'])

if __name__ == '__main__':
    make_data()