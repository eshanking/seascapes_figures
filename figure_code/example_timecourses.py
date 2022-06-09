import os
import matplotlib.pyplot as plt
import numpy as np
import pickle
from seascapes_figures.utils import results_manager, plotter

def load_exp_info(info_path):
    return pickle.load(open(info_path, 'rb'))

# def make_figure(roc_exp,adh_exp):
def make_figure(roc_info_path,adh_info_path):
    
    fig,ax = plt.subplots(nrows=3,ncols=4,figsize=(6,4))
    labelsize = 10
    #%% ROC data
    # suffix = '11012021_0000' # lab machine
    # # suffix = '10272021_0001' # macbook
    # data_folder = 'results_' + suffix
    # exp_info_file = 'experiment_info_' + suffix + '.p'
    
    roc_exp = load_exp_info(roc_info_path)
    adh_exp = load_exp_info(adh_info_path)

    # data_folder = roc_exp.results_path
    # exp_info_file = roc_exp.experiment_info_path
    
    exp_folders,exp_info = results_manager.get_experiment_results(exp=roc_exp)
    # max_cells = exp_info.populations[0].max_cells
    # n_sims = exp_info.n_sims
    # k_abs = exp_info.slopes
    
    exp_num = 1
    pop_roc = exp_info.populations[0]
    

    
    # find one extinct and one survived simulation
    
    # k = 0
    found_extinct = False
    found_survived = False

    n_sims = roc_exp.n_sims
    
    while found_extinct == False or found_survived == False:
        # while exp_num < len(exp_folders):
        exp = exp_folders[exp_num]
        sim_files = os.listdir(path=exp)
        sim_files = sorted(sim_files)
        k = 0

        while k < n_sims:
            sim = sim_files[k]
            sim = exp + os.sep + sim
            data = results_manager.get_data(sim)
            # data = data[:,0:-1]
            # data_t = data[-1,:]
            data_t = data['counts']
            data_t = data_t[-1,:]
            if any(data_t >= 1):
                num_survived = k
                found_survived = True
            else:
                num_extinct = k
                # exp_extinct = exp_num
                found_extinct = True
            k+=1
            # exp_num += 1


    #%%
    # plot ROC survived timecourse
    
    sim = sim_files[num_survived]
    sim = exp + os.sep + sim
    data = results_manager.get_data(sim)
    dc = data['drug_curve']
    # data = data[:,0:-1]
    # data = data/np.max(data)
    data_t = data['counts']
    tcax = ax[0,0]
    drug_kwargs = {'alpha':1,
                   'color':'black',
                   'linewidth':1,
                   'linestyle':'-'
                   # 'label':'Drug Concentration ($\u03BC$M)'
                   }
    label_kwargs = {'align':False}
    select_labels = [7,15]
    label_xpos = [900,2500]
    
    data_t = data_t/np.max(data_t)

    p = roc_exp.populations[0]
    
    tcax,t = p.plot_timecourse_to_axes(data_t,
                                        tcax,
                                        # drug_curve=dc,
                                        drug_ax_sci_notation=True,
                                        drug_kwargs=drug_kwargs,
                                        legend_labels=True,
                                        linewidth=1.5,
                                        drug_curve_label = '',
                                        labelsize = labelsize,
                                        label_lines = True,
                                        select_labels = select_labels,
                                        label_xpos=label_xpos,
                                        label_kwargs=label_kwargs
                                        )
    drug_ax = ax[1,0]
    drug_ax.plot(dc,color='black',linewidth=1)
    # drug_ax.set_yticklabels('')
    # drug_ax.set_yticks([])
    
    #%% plot ROC extinct timecourse
    
    sim = sim_files[num_extinct]
    sim = exp + os.sep + sim
    data = results_manager.get_data(sim)
    dc = data['drug_curve']
    # data = data[:,0:-1]
    # data = data/np.max(data)
    data_t = data['counts']
    tcax = ax[0,1]
    drug_kwargs = {'alpha':1,
                   'color':'black',
                   'linewidth':1
                   # 'label':'Drug Concentration ($\u03BC$M)'
                   }
    
    data_t = data_t/np.max(data_t)

    select_labels = [2]
    label_xpos = [500]
    
    label_kwargs = {'align':True}

    tcax,t = p.plot_timecourse_to_axes(data_t,
                                        tcax,
                                        # drug_curve=dc,
                                        drug_ax_sci_notation=True,
                                        drug_kwargs=drug_kwargs,
                                        legend_labels=True,
                                        linewidth=1.5,
                                        labelsize = labelsize,
                                        label_lines = True,
                                        select_labels = select_labels,
                                        label_xpos=label_xpos,
                                        label_kwargs=label_kwargs
                                        )
    drug_ax = ax[1,1]
    drug_ax.plot(dc,color='black',linewidth=1)
    #%% nonadherance data
    
    # suffix = '11042021_0000' # lab machine
    # exp_num = 2 # lab machine
    # # exp_num = 0 # macbook
    # # suffix = '11082021_0000' # macbook
    # data_folder = 'results_' + suffix
    # exp_info_file = 'experiment_info_' + suffix + '.p'
    
    data_folder = adh_exp.results_path
    exp_info_file = adh_exp.experiment_info_path
    
    exp_folders,exp_info = results_manager.get_experiment_results(exp=adh_exp)
    max_cells = exp_info.populations[0].max_cells
    n_sims = exp_info.n_sims
    
    pop = exp_info.populations[0]
    pop_adh = pop
    n_timestep = pop.n_timestep
    
    exp_num = 3
    # exp = exp_folders[exp_num]
    
    sim_files = os.listdir(path=exp)
    sim_files = sorted(sim_files)
    
    # find one extinct and one survived simulation
    
    k = 0
    found_extinct = False
    found_survived = False
    
    # while found_extinct == False or found_survived == False:
    #     sim = sim_files[k]
    #     sim = exp + os.sep + sim
    #     data = results_manager.get_data(sim)
    #     # data = data[:,0:-2]
    #     data_t = data['counts']
    #     data_t = data_t[-1,:]
    #     if any(data_t >= 1):
    #         num_survived = k
    #         found_survived = True
    #     else:
    #         num_extinct = k
    #         found_extinct = True
    #     k+=1

    k = 0
    while (found_extinct == False or found_survived == False) and k < n_sims:
        # while exp_num < len(exp_folders):
        exp = exp_folders[exp_num]
        sim_files = os.listdir(path=exp)
        sim_files = sorted(sim_files)

        sim = sim_files[k]
        sim = exp + os.sep + sim
        data = results_manager.get_data(sim)
        # data = data[:,0:-1]
        # data_t = data[-1,:]
        data_t = data['counts']
        data_t = data_t[-1,:]
        if any(data_t >= 1):
            num_survived = k
            found_survived = True
        else:
            num_extinct = k
            # exp_extinct = exp_num
            found_extinct = True
        k+=1
    
    if found_extinct == False or found_survived == False:
        print(found_extinct)
        print(found_survived)
        raise ValueError('No suitable example traces found')
    #%% plot nonadherance survived timecourse
    
    sim = sim_files[num_survived]
    sim = exp + os.sep + sim
    data = results_manager.get_data(sim)
    dc = data['drug_curve']
    survived_schedule = data['regimen']
    
    tcax = ax[0,2]
    drug_kwargs = {'alpha':1,
                   'color':'black',
                   'linewidth':1
                   # 'label':'Drug Concentration ($\u03BC$M)'
                   }
    label_kwargs = {'align':True}
    tc = data['counts']
    tc = tc/np.max(tc)
    
    select_labels = [0,2]
    label_xpos = [150,550]
    
    tcax,drug_ax = p.plot_timecourse_to_axes(tc,
                                            tcax,
                                            # drug_curve=dc,
                                            drug_ax_sci_notation=True,
                                            drug_kwargs=drug_kwargs,
                                            legend_labels=True,
                                            legend_size=16,
                                            linewidth=1.5,
                                            drug_curve_label = '',
                                            labelsize = labelsize,
                                            label_lines=True,
                                            select_labels=select_labels,
                                            label_xpos=label_xpos,
                                            label_kwargs=label_kwargs
                                            )
    drug_ax = ax[1,2]
    drug_ax.plot(dc,color='black',linewidth=1)
    
    #%% plot nonadherance extinct timecourse
    
    sim = sim_files[num_extinct]
    sim = exp + os.sep + sim
    data = results_manager.get_data(sim)
    dc = data['drug_curve']
    extinct_schedule = data['regimen']
    
    tcax = ax[0,3]
    drug_kwargs = {'alpha':1,
                   'color':'black',
                   'linewidth':1
                   # 'label':'Drug Concentration ($\u03BC$M)'
                   }
    
    tc = data['counts']
    tc = tc/np.max(tc)
    
    select_labels = [0]
    label_xpos = [100]
    
    tcax,drug_ax = p.plot_timecourse_to_axes(tc,
                                                tcax,
                                                # drug_curve=dc,
                                                drug_ax_sci_notation=True,
                                                drug_kwargs=drug_kwargs,
                                                legend_labels=True,
                                                linewidth=1.5,
                                                labelsize = labelsize,
                                                legend_size=16,
                                                select_labels=select_labels,
                                                label_xpos=label_xpos,
                                                label_lines=True
                                                )
    drug_ax = ax[1,3]
    drug_ax.plot(dc,color='black',linewidth=1)
    #%%  plot dose times
    
    x = np.arange(n_timestep)
    
    timescale = pop.dose_schedule/pop.timestep_scale
    
    ax[2,2].plot(x,survived_schedule,linewidth=0.2,color='black')
    ax[2,3].plot(x,extinct_schedule,linewidth=0.2,color='black')
    
    #%% Adjust positions
    
    for col in range(0,4):
        ax[0,col] = p.shrinky(ax[0,col],0.04)
        ax[1,col] = p.shrinky(ax[1,col],0.11)
        # ax[1,col] = p.shifty(ax[1,col],-0.04)
        ax[1,col].ticklabel_format(style='scientific',axis='y',
                                              scilimits=(0,3))
        
        ax[2,col] = p.shrinky(ax[2,col],0.2)
        ax[1,col] = p.shifty(ax[1,col],-0.02)
        ax[2,col] = p.shifty(ax[2,col],0.10)
    
    # shift right column to the right to make room for axis labels
    for col in range(2,4):
        for row in range(0,3):
            ax[row,col] = p.shiftx(ax[row,col],0.1)
        
            
    ax[2,2].spines['top'].set_visible(False)
    ax[2,3].spines['top'].set_visible(False)
    ax[2,2].spines['left'].set_visible(False)
    ax[2,3].spines['left'].set_visible(False)
    ax[2,2].spines['right'].set_visible(False)
    ax[2,3].spines['right'].set_visible(False)
    
    for row in range(0,3):
        ax[row,1].set_yticks([])
        ax[row,3].set_yticks([])
    
    ax[2,2].set_yticks([])
    
    ax[2,0].remove()
    ax[2,1].remove()
    
    #%% Adjust x axis tick labels
    
    # ax[0,0] = pop_roc.x_ticks_to_days(ax[0,0])
    time_scale = 24/pop_roc.timestep_scale
    ax[0,0].set_xticks([0,250*time_scale,500*time_scale])
    ax[0,0].set_xticklabels([0,250,500])

    # xlim = [0,500*time_scale]

    # ax[0,0].set_xlim(xlim)
    xt = ax[0,0].get_xticks()    
    xl = ax[0,0].get_xticklabels()
    xlim = ax[0,0].get_xlim()

    ax[0,1].set_xticks(xt)
    ax[0,1].set_xticklabels(xl)
    ax[0,1].set_xlim(xlim)
    
    for col in range(0,2):
        ax[1,col] = pop_roc.x_ticks_to_days(ax[1,col])
        ax[1,col].set_xticks(xt)
        ax[1,col].set_xticklabels(xl)
        ax[1,col].set_xlim(xlim)
        ax[1,col].set_xlabel('Days')
    

    for col in range(2,4):
        for row in range(0,3):
            ax[row,col] = pop_adh.x_ticks_to_days(ax[row,col])

    xt = ax[1,2].get_xticks()
    xl = ax[1,2].get_xticklabels()
    xlim = ax[1,2].get_xlim()

    ax[0,2].set_xticks(xt)
    ax[0,2].set_xticklabels(xl)
    ax[0,2].set_xlim(xlim)

    ax[0,3].set_xticks(xt)
    ax[0,3].set_xticklabels(xl)
    ax[0,3].set_xlim(xlim)

    ax[1,0].ticklabel_format(style='scientific',axis='y',scilimits=(0,1))
    ax[1,2].ticklabel_format(style='scientific',axis='y',scilimits=(0,1))
    
    #%% Add axis labels
    
    ax[0,0].set_ylabel('Proportion',fontsize=8,labelpad=0.8)
    ax[1,0].set_ylabel('Concentration ($\u03BC$M)',fontsize=8)
    
    ax[0,2].set_ylabel('Proportion',fontsize=8,labelpad=0.8)
    ax[1,2].set_ylabel('Concentration ($\u03BC$M)',fontsize=8)
    
    #%% add panel labels
    
    alphabet = ['a','b','','d','e','f','','h','i','j']
    
    k = 0
    
    for row in range(2):
        for col in range(2):
            ax[row,col].annotate(alphabet[k],(0,1.07),xycoords='axes fraction')
            k+=1
            
    for row in range(0,2):
        for col in range(2,4):
            ax[row,col].annotate(alphabet[k],(0,1.1),xycoords='axes fraction')
            k+=1
    
    ax[2,2].annotate(alphabet[k],(0,1.4),xycoords='axes fraction')
    k+=1
    ax[2,3].annotate(alphabet[k],(0,1.4),xycoords='axes fraction')

    ax[1,0].annotate('c',(0,1.47),xycoords='axes fraction')
    ax[1,2].annotate('g',(0,1.47),xycoords='axes fraction')
    
    # fig.tight_layout()

    results_manager.save_fig(fig,'example_timecourses.pdf',bbox_inches='tight')
    
# if __name__ == '__main__':
#     make_figure()     

make_figure(rip,aip)                       