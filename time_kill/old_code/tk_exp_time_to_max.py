# %%
from fears.utils import AutoRate
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy

mic = 0.08 # mic of strain 0

# sample_times = [0,30,60,90,120,150,180,220,240,270,300,330]
sample_times = []

# sample_times = np.delete(sample_times,4)

exp_folder = 'tk_02282023'

exp_files = os.listdir(exp_folder)

exp_files = [f for f in exp_files if f[-5:] == '.xlsx']

exp_files.sort()

# exp_files = exp_files

n_samples = len(exp_files)

# Load and normalize each data file

plate_list = []

row_list = ['A','B','C','D','E','F','G','H']

# drug_conc = [0,0.25,0.5,1,2,4,8]
drug_conc = [0,'$10^{-2}$','$10^{-1}$','1','10','$10^{2}$','$10^{3}$']
drug_conc = [str(dc) for dc in drug_conc]

# Helper functions

def rolling_average(x,N):
    
    indx = 0
    res = []
    while indx+N < len(x):
        xt = x[indx:indx+N]
        res.append(np.nanmean(xt))
        indx+=1
    
    for i in range(len(x[indx:])):
        res.append(np.nanmean(x[indx+i:]))

    res = np.array(res)
    x = np.array(x)
    res[x == 0] = 0

    return res

"""
# Model adapted from
Li RC. Simultaneous pharmacodynamic analysis of the lag and bactericidal phases
exhibited by beta-lactams against Escherichia coli. Antimicrob Agents 
Chemother. 1996;40(10):2306-2310. doi:10.1128/AAC.40.10.2306
"""

# def time_kill_model(t,kss,alpha,y0):
#     # Growth rate constant of control arbitrarily set to zero
#     res = []
#     for tt in t:
#         res.append(y0 - kss*tt + (kss/alpha)*(1-np.exp(-1*alpha*tt)))
#     return res

def time_kill_model(t,alpha,A,l): # decreasing signal
    res = []
    for tt in t:
        res.append(A*np.exp(alpha*(tt-l)))
    return res

def est_time_kill(xdata,ydata,debug=False):


    # alpha,A,l

    A_est = ydata[0]
    alpha_est = 0
    p0 = [alpha_est,A_est,0]
    # bounds = [[-np.inf,A_est-0.2*(np.abs(A_est)),-np.inf],
    #     [np.inf,A_est+0.2*(np.abs(A_est)),np.inf]]

    popt,pcov = scipy.optimize.curve_fit(time_kill_model,xdata,ydata,p0=p0)

    if debug:

        y_opt = time_kill_model(xdata,popt[0],popt[1],popt[2])
        
        fig,ax = plt.subplots()
        ax.plot(xdata,y_opt)
        ax.plot(xdata,ydata)
        ax.set_title(popt[0])
        # ax.set_ylim(-1000,50000)

    return popt

def pharmacodynamic_curve(c, gmax, gmin, mic, k):
    """pharmacodynamic model adapted from Foerster et al.

    Foerster, S., Unemo, M., Hathaway, L.J. et al. Time-kill curve analysis and
    pharmacodynamic modelling for in vitro evaluation of antimicrobials against Neisseria
    gonorrhoeae . BMC Microbiol 16, 216 (2016). 
    https://doi.org/10.1186/s12866-016-0838-9

    Args:
        c (float): drug concentration
        gmax (float): max growth rate
        gmin (float): min growth rate
        mic (float): estimated minimum inhibitory concentration
        k (float): hill coefficient
    """
    g = []
    for c_t in c:
        g.append(gmax - (((gmax-gmin)*(c_t/mic)**k)/((c_t/mic)**k-(gmin/gmax))))
    
    return g

def est_pharm_curve(xdata,ydata,debug=False):


    # gmax,gmin,mic,k

    g_max_est = ydata[0]
    gmin_est = np.min(ydata)

    p0 = [g_max_est,gmin_est,1,0.1]

    # p0 = [1,g_drugless_est,-1,0]
    # bounds = ([-np.inf,g_drugless_est-g_drugless_est*0.05,-2,0],
    #           [np.inf,g_drugless_est+g_drugless_est*0.05,-0.1,np.inf])

    popt,pcov = scipy.optimize.curve_fit(pharmacodynamic_curve,xdata,ydata,p0=p0)

    if debug:
        
        dc_fit = np.logspace(-3,3)
        y_opt = pharmacodynamic_curve(dc_fit,popt[0],popt[1],popt[2],popt[3])
        
        fig,ax = plt.subplots()
        ax.plot(dc_fit,y_opt)
        ax.plot(xdata,ydata,'o')

        ax.set_xscale('log')

    return popt

#%% Background normalization and data cleaning

od_data_list = []
cur_col = 1

for ef in exp_files:

    path_t = os.getcwd() + os.sep + exp_folder + os.sep + ef
    p = AutoRate.Plate(path_t)

    od_data = p.parse_data_file(path_t,data_start='OD600')
    od_data_list.append(od_data)
    
    # normalize to bottom row
    # estimate background and exclude 6,7,9,10,11,12 because of contamination
    # bg_key = 'H1'
    # bg_data = np.array(p.data[bg_key])

    # for col in [2,3,4,5,8]:
    #     bg_key = row_list[-1] + str(col)
    #     bg_data += np.array(p.data[bg_key])
    
    # bg_data = bg_data/5

    p.data = p.data.drop(p.data.index[0:2]) # remove first datapoint because of a systematic drop in fluorescence

    for col in range(1,cur_col+1):
        
        bg_key = row_list[-1] + str(col)
        bg_data = p.data[bg_key]

        for row in row_list[0:-1]:

            key = row + str(col)
            ts = p.data[key]

            bg_data[bg_data==0] = 1

            ts[ts=='OVER'] = np.nan
            ts_norm = np.divide(ts,bg_data)

            p.data[key] = ts_norm
    
    for col in range(cur_col+1,13):
        
        # all fluorescence in columns after the column corresponding to the current
        # column index are set to zero
        for row in row_list[0:-1]:
            key = row + str(col)

            p.data[key] = p.data[key]*0

    # add start time to time
    dt = p.get_start_time()
    start_time = 60*((60*dt.hour) + dt.minute) + dt.second
    sample_times.append(start_time/60)
    # print(start_time)
    p.data['Time [s]'] = p.data['Time [s]'] + start_time

    od_data['Time [s]'] = od_data['Time [s]'] + start_time

    plate_list.append(p)
    cur_col+=1

# %% Combine and sort timecourses

data = plate_list[0].data
od_data = od_data_list[0]

for od_d in od_data_list[1:]:
    od_data = pd.concat((od_data,od_d))

for p in plate_list[1:]:
    data = pd.concat((data,p.data))

data = data.sort_values(by='Time [s]')

data['Time [s]'] = data['Time [s]'] - np.min(data['Time [s]'])

od_data['Time [s]'] = od_data['Time [s]'] - np.min(od_data['Time [s]'])

sample_times = np.array(sample_times) - np.min(sample_times)

# Estimate peak fluorescence for each drug concentration at each sample time

#%% Plot raw fluorescence timecourse
fig,ax_list = plt.subplots(nrows=4,ncols=2,figsize=(8,10))

ax_list_t = ax_list.reshape(-1)

ax_indx = 0

time = np.array(data['Time [s]'])

cmap = mpl.colormaps['viridis']

for row in row_list[0:-1]:
    col_indx = 0
    ax = ax_list_t[ax_indx]
    # for col in range(1,13):
    for col in [1,2,3,4,5,8]:

        key = row + str(col)
        
        if row == 'A':
            label=round(sample_times[col-1])
        else:
            label=None

        sample_time = sample_times[col-1]*60 - 30*60
        time_t = time - sample_time
        # time_t = time

        ts = np.array(data[key]).astype(float)

        # ts = rolling_average(ts,10)

        # data[key] = ts

        ax.plot(time_t,ts,label=label,color=cmap(col/6),linewidth=2)

        # indx = np.argwhere(time_t>=sample_time)[0][0]

        # ax.plot(time_t[indx],ts[indx],'*',color='black')

        # max_fluor = np.nanmax(ts)
        # indx = np.argwhere(ts==max_fluor)[0][0]

        # ax.plot(time_t[indx],ts[indx],'*',markersize=5,color='black')
        

        # ax.set_title(str(sample_times[ax_indx]))

        # ax.set_xlim(0,17000)
        col_indx+=1
    ax.set_xlim(-100,80000)
    ax.set_title(drug_conc[ax_indx])
    ax_indx+=1

for row in range(4):
    for col in range(2):
        pos = ax_list[row,col].get_position()
        pos.y0 = pos.y0 - 0.02*row
        pos.y1 = pos.y1 - 0.02*row
        ax_list[row,col].set_position(pos)

fig.legend(ncol=6,frameon=False,loc=(0.1,0))

fig.savefig('all_fluorescence_timecouses.pdf',bbox_inches='tight')
# %%

# Get time to peak
time_dict = {}

time_vect = np.array(data['Time [s]'])

for row in row_list[0:-2]: # for each drug conc

    time_t = []
    for col in range(1,13):
        sample_time = sample_times[col-1]*60
        time_vect_t = time_vect - sample_time
        key = row+str(col)
        ts = np.array(data[key])
        # ts = rolling_average(ts,10)

        if row == 'G':
            # max_fluor = np.nanmax(ts[200:])
            ts = ts[200:]
        max_fluor = np.nanmax(ts)
        indx = np.argwhere(ts==max_fluor)[0][0]

        if row == 'G':
            indx += 200

        time_t.append(time_vect_t[indx])
    
    time_dict[row] = np.array(time_t)   

# plot raw results

fig2,ax = plt.subplots()

row_indx = 0
for row in row_list[0:-2]:
    line = ax.plot(sample_times[0:n_samples],time_dict[row][0:n_samples],label=drug_conc[row_indx])
    c = line[0].get_color()
    ax.plot(sample_times[0:n_samples],time_dict[row][0:n_samples],'*',color=c)
    row_indx+=1

# ax.set_xlim(150,300)

ax.legend(frameon=False,title='xMIC')

ax.set_ylabel('Time to max',fontsize=14)
ax.set_xlabel('Sample time (min)',fontsize=14)

fig2.savefig('max_fluorescence_over_time.pdf',bbox_inches='tight')

#%% Fit exponentials

alpha = []

for row in row_list[0:-2]:

    # norm_const = np.min(diff_dict[row])

    ydata = time_dict[row]
    ydata = ydata - np.min(ydata)
    # ydata = ydata/np.max(ydata)

    # if ydata[-1] > ydata[0]:
    #     ydata_t = 1-ydata
    # else:
    #     ydata_t = ydata

    # if row == 'F':
    ydata = ydata[0:-2]
    sample_times_t = sample_times[0:-2]
    # else:
    #     sample_times_t = sample_times

    popt = est_time_kill(sample_times_t,ydata,debug=True)

    # if ydata[-1] > ydata[0]:
    #     alpha.append(-popt[0])
    # else:
    alpha.append(-popt[0])

#%% plot dose-response curve

dc = [-3,-2,-1,0,1,2]

fig3,ax = plt.subplots()

ax.plot(dc,alpha,'o')

# %% fit logistic curve

dc_t = [10**c for c in dc]

popt = est_pharm_curve(dc_t,alpha,debug=True)

#%% scale dose-response curve

fig,ax = plt.subplots()

dc_fit = np.logspace(-3,3)

y_opt = pharmacodynamic_curve(dc_fit,popt[0],popt[1],popt[2],0.1)

# indx = np.argwhere(dc_fit>0.08)[0][0]
# y_opt = y_opt - y_opt[indx]

ax.plot(dc_fit,y_opt)
ax.set_xscale('log')
ax.scatter(dc_t,alpha,color='black',marker='x')

ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
ax.tick_params(axis='both', labelsize=12)

ax.set_ylabel('Growth rate (AU)',fontsize=14)
ax.set_xlabel('Drug concentration (ug/mL)',fontsize=14)
ax.set_title('Genotype 0 time-kill assay',fontsize=14)

fig.savefig('time-kill.pdf',bbox_inches='tight')

# %%
