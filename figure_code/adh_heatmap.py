import os
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoLocator
import numpy as np
from seascapes_figures.utils import results_manager, plotter, dir_manager
import pandas as pd
import pickle

def make_fig(exp=None,exp_info_path=None):

    if exp is None:
        exp = pickle.load(open(exp_info_path,'rb'))

    exp_folders,exp_info = results_manager.get_experiment_results(exp=exp)

    n_sims = exp_info.n_sims
    p_drop = exp_info.prob_drops

    exp_folders.reverse()
    p_drop = np.flip(p_drop)

    pop = exp_info.populations[0]
    
    exp = exp_folders[2]

    p_drop_t = exp[exp.find('=')+1:]
    p_drop_t = p_drop_t.replace(',','.')
    p_drop_t = float(p_drop_t)
    
    num = np.argwhere(p_drop == p_drop_t)
    num = num[0,0]
    
    sim_files = os.listdir(path=exp)
    sim_files = sorted(sim_files)

    pop = exp_info.populations[0]

    gap = int(pop.dose_schedule/pop.timestep_scale)
    n_scheduled_doses = int(np.ceil(pop.n_timestep/gap))

    while k < len(sim_files):
        sim = sim_files[k]
        sim = exp+ os.sep + sim
        data_dict = results_manager.get_data(sim)
        counts = data_dict['counts']
