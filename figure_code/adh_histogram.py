import os
import matplotlib.pyplot as plt
# from matplotlib.ticker import AutoLocator
import matplotlib as mpl
import numpy as np
from fears.utils import results_manager
import pickle
from scipy import stats
import seaborn as sns
import pandas as pd

def make_fig(exp=None,exp_info_path=None):

    if exp is None:
        exp = pickle.load(open(exp_info_path,'rb'))

    exp_folders,exp_info = results_manager.get_experiment_results(exp=exp)

    n_sims = exp_info.n_sims
    p_drop = exp_info.prob_drops

    exp_folders.reverse()
    p_drop = np.flip(p_drop)

    pop = exp_info.populations[0]
    
    # exp = exp_folders[2]

    pop = exp_info.populations[0]

    gap = int(pop.dose_schedule/pop.timestep_scale)     
    if pop.dwell:
        dwell_doses = int(np.floor(pop.dwell_time/gap))
    else:
        dwell_doses = 0

    n_scheduled_doses = int(np.floor(pop.n_timestep/gap))-dwell_doses

    max_doses = int(np.sum(pop.impulses))

    n_subplots = len(exp_folders)-1
    fig,ax_list = plt.subplots(ncols=n_subplots,nrows=2,figsize=(10,5))

    # fig2,ax_tb = plt.subplots(ncols=n_subplots,figsize=(10,3))

    # k=0

    # for exp in [exp_folders[2]]:
    num = 0
    for exp in exp_folders[:-1]:
        ax = ax_list[0,num]
        extinct_sched = np.zeros(n_scheduled_doses)
        extinct_sched = np.array([extinct_sched])
        survived_sched = np.zeros(n_scheduled_doses)
        survived_sched = np.array([survived_sched])
        sim_files = os.listdir(path=exp)
        sim_files = sorted(sim_files)

        p_drop_t = exp[exp.find('=')+1:]
        p_drop_t = p_drop_t.replace(',','.')
        p_drop_t = float(p_drop_t)
        
        num = np.argwhere(p_drop == p_drop_t)
        num = num[0,0]
        k=0
        while k < len(sim_files):

            sim = sim_files[k]
            sim = exp+ os.sep + sim
            data_dict = results_manager.get_data(sim)
            counts = data_dict['counts']
            counts = np.sum(counts,axis=1)
            # dc = data_dict['drug_curve']
            regimen = data_dict['regimen']
            dose_schedule = compute_doses(regimen,pop)
            dose_schedule = dose_schedule[dwell_doses:]
            dose_schedule = np.array([dose_schedule])
            # dose_schedule = dose_schedule[dwell_doses:]

            if counts[-1]<1:
                extinct_sched = np.concatenate((extinct_sched,dose_schedule),axis=0)
            else:
                survived_sched = np.concatenate((survived_sched,dose_schedule),axis=0)
            k+=1

        extinct_sched = extinct_sched[1:,:]
        survived_sched = survived_sched[1:,:]

        # fig,ax = plt.subplots(2,1,figsize=(6.25,7.75),sharex=True)
        # fig,ax = plt.subplots(figsize=(4,3))
        cmap = mpl.colors.ListedColormap(['cornflowerblue','w'])

        # aspect = n_scheduled_doses/n_sims
        # aspect_surv = n_scheduled_doses/survived_sched.shape[0]

        # ax[0].imshow(extinct_sched,cmap=cmap)
        # ax[1].imshow(survived_sched,cmap=cmap)

        n_extinct = extinct_sched.shape[0]
        n_survived = survived_sched.shape[0]

        survived_hist = np.sum(survived_sched,axis=0)
        sh_counts = survived_hist
        survived_hist = survived_hist/n_survived
        extinct_hist = np.sum(extinct_sched,axis=0)
        ext_counts = extinct_hist
        extinct_hist = extinct_hist/n_extinct

        dose_num = np.arange(len(survived_hist)) + 1
        # p = (1-float(p_drop_t))*np.ones(len(dose_num))

        ratio = np.divide(extinct_hist,survived_hist)

        ax.bar(dose_num,ratio,width=1,color='red',alpha=0.6,edgecolor='w',
            align='edge')
        ax.plot(dose_num,np.ones(len(dose_num)),'--',color='black')
    
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)

        
        ax.set_xlabel('Scheduled dose',fontsize=12)

        ax.tick_params(axis='x', labelsize=12)
        ax.tick_params(axis='y', labelsize=12)

        ax.set_xlim(1,max_doses+1)
        ax.set_ylim(0,3.5)

        ax.set_title('$p_{forget} = $' + str(p_drop_t))

        ax.text(x=3,y=0.3,s='$n_{success} = $' + str(n_extinct),
                backgroundcolor='gray',color='w')

        # Statistical analysis

        # res = []

        # sig_level = 0.05/max_doses # bonferroni correction

        # for i in range(len(ext_counts)):
        #     ext_taken = ext_counts[i]
        #     ext_not_taken = n_extinct-ext_taken

        #     surv_taken = sh_counts[i]
        #     surv_not_taken = n_survived-surv_taken

        #     if surv_taken + ext_taken > 0:
        #         obs = np.array([[ext_taken,ext_not_taken],[surv_taken,surv_not_taken]])
                
        #         zeros = np.argwhere(obs==0)

        #         if not (zeros.shape[0]>0):

        #             p = stats.chi2_contingency(obs)[1]

        #             res.append(p)

        #             if p < sig_level:
        #                 ax.text(dose_num[i]+0.1,ratio[i],'*',fontsize='15')

        # Time between first and second dose

        ax = ax_list[1,num]

        time_between_surv = []
        time_between_ext = []

        for s in survived_sched:
            # print(s)
            doses = np.argwhere(s==1)
            if len(doses) > 1:
                time_between_surv.append((doses[1]-doses[0])[0])
        for e in extinct_sched:
            doses = np.argwhere(e==1)
            if len(doses) > 1:
                time_between_ext.append((doses[1]-doses[0])[0])

        d = {'survived':time_between_surv,
             'extinct':time_between_ext}

        # df = pd.DataFrame(d)

        sns.stripplot(ax=ax,data=[time_between_surv,time_between_ext],size=3,
                      jitter=0.2)
        
        # ax.tight_layout()
        ax.set_xticklabels(['failure','success'])

        ax.set_yticks([0,5,10,15,20])
        ax.set_ylim(0,20)

        ax.tick_params(axis='x', labelsize=12)
        ax.tick_params(axis='y', labelsize=12)

        # ax.set_ylabel('Time between first doses',fontsize=12)

        num+=1

    ax_list[0,0].set_ylabel('Odds (success/failure)',fontsize=12)
    ax_list[1,0].set_ylabel('Time between first doses',fontsize=12)
    # fig2.tight_layout()
    fig.tight_layout()
    fig.savefig('figures/dose_ratio_barplot.pdf',bbox_inches='tight')
    # return survived_sched,extinct_sched
    # return df



def compute_doses(regimen,pop):

    dose_schedule = pop.dose_schedule/pop.timestep_scale
    n_doses = int(np.floor(pop.n_timestep/dose_schedule))

    doses = np.zeros(n_doses)

    dose_num = 0

    for i in range(n_doses):
        if regimen[int(i*dose_schedule)] == 1:
            doses[dose_num] = 1
        dose_num += 1
    
    return doses
    
