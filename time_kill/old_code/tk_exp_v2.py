#%%
import pickle
from fears.utils import AutoRate
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.signal import butter, lfilter, freqz
import scipy

data_file = 'EK_AB_20230321_095345.xlsx'

row_list = ['A','B','C','D','E','F','G','H']
col_list = [2,3,4,5,6,7,8,9,10,11]
col_list = [str(col) for col in col_list]

sample_times = [0,30,90,210,390]
# 12:49:21
st_data_cleaning = [32*60+33,92*60+27,215*60+25,392*60+14]

dc = [0,10**-2,10**-1,10**0,10,10**2,10**3,0]
mic = 0.1

dc = [c*mic for c in dc]

# load normalizer object

filehandler = open('normalizer.p', 'rb') 
norm = pickle.load(filehandler)

#%% Stitcher function to stich multiple scans together

def stitcher(exp_folder,background=None,bg_row='H'):

    plate_list = []
    od_data_list = []

    exp_folder = 'tk_02282023'

    exp_files = os.listdir(exp_folder)

    exp_files = [f for f in exp_files if f[-5:] == '.xlsx']

    cur_col = 1

    for ef in exp_files:
        path_t = os.getcwd() + os.sep + exp_folder + os.sep + ef
        p = AutoRate.Plate(path_t)

        od_data = p.parse_data_file(path_t,data_start='OD600')
        od_data_list.append(od_data)
        
        p.data = p.data.drop(p.data.index[0:2])


        for col in range(1,cur_col+1):
            
            bg_key = bg_row + str(col)
            bg_data = p.data[bg_key]

            for row in row_list[0:-1]:

                key = row + str(col)
                ts = p.data[key]

                bg_data[bg_data==0] = 1

                ts[ts=='OVER'] = np.nan
                if background == 'subtract':
                    ts = ts - bg_data
                elif background == 'divide':
                    ts = np.divide(ts,bg_data)

                p.data[key] = ts
        dt = p.get_start_time()
        start_time = 60*((60*dt.hour) + dt.minute) + dt.second
        sample_times.append(start_time/60)
        # print(start_time)
        p.data['Time [s]'] = p.data['Time [s]'] + start_time

        od_data['Time [s]'] = od_data['Time [s]'] + start_time

        plate_list.append(p)
        cur_col+=1

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
    data['est. sample times']
                    



#%% helper functions
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

def est_time_kill(xdata,ydata,debug=False):
    # alpha,A,l
    
    indx = np.argwhere(np.abs(ydata)==np.max(np.abs(ydata)))[0][0]
    A_est = ydata[indx]

    alpha_est = 0
    p0 = [alpha_est,A_est,0]
    # bounds = [[-np.inf,A_est-0.2*(np.abs(A_est)),-np.inf],
    #     [np.inf,A_est+0.2*(np.abs(A_est)),np.inf]]

    popt,pcov = scipy.optimize.curve_fit(time_kill_model,xdata,ydata,p0=p0)

    if debug:

        y_opt = time_kill_model(xdata,popt[0],popt[1],popt[2])
        
        fig,ax = plt.subplots()
        ax.plot(xdata,y_opt,label='est')
        ax.scatter(xdata,ydata,marker="*",label='data')
        ax.set_title(popt[0])
        ax.legend()
        # ax.set_ylim(-1000,50000)

    return popt

# def time_kill_model(t,alpha,A,l): # decreasing signal
#     res = []
#     for tt in t:
#         res.append(A*np.exp(alpha*(tt-l)))
#     return res

# def time_kill_model(t,alpha,A,l,y0):
#     res = []
#     for tt in t:
#         res.append(y0 + A/(1+np.exp(alpha*(tt-l))))
#     return res

def time_kill_model(t,alpha,A,l):
    res = []
    for tt in t:
        res.append(A*(1-np.exp(-alpha*(tt-l)/A)))
    return res

def butter_lowpass(cutoff, fs, order=5):
    return butter(order, cutoff, fs=fs, btype='low', analog=False)

def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y

order = 1
fs = 0.01      # sample rate, Hz
cutoff = 0.0001  # desired cutoff frequency of the filter, Hz

# %% data cleaning

p = AutoRate.Plate(data_file)
time_vect = np.array(p.data['Time [s]'])

od_data = p.parse_data_file(data_file,data_start='OD600')

data_cleaning_indx = []
for tt in st_data_cleaning:
    indx = np.argwhere(time_vect>=tt)[0][0]
    p.data = p.data.drop(p.data.index[indx-1:indx+2])
    data_cleaning_indx.append(indx)
    time_vect = np.array(p.data['Time [s]'])

p.data = p.data.drop(p.data.index[-1]) # remove last datapoint because of incomplete data
time_vect = np.array(p.data['Time [s]'])

fig,ax_list = plt.subplots(nrows=4,ncols=2,figsize=(8,10))
ax_list_t = ax_list.reshape(-1)
cmap = mpl.colormaps['viridis']

norm_data = {}
od_data = p.g

row_indx = 0
for row in row_list:
    ax = ax_list_t[row_indx]
    col_indx = 0
    max_fluor = 0
    for col in col_list:
        key = row+col
        ts = np.array(p.data[key])
        conc = dc[row_indx]
        if conc>0:
            conc = np.log10(conc)
            for indx in range(len(ts)):
                if time_vect[indx] >= np.max(np.array(norm.data['Time [s]'])):
                    time_t = np.max(np.array(norm.data['Time [s]']))
                else:
                    time_t = time_vect[indx]
                
                if time_t > 11000 and time_t < 16000:
                    time_t = 16000
                norm_factor = norm.get_ctx_norm_factor(conc,time_t)
                ts[indx] = ts[indx]/norm_factor
        # ts = butter_lowpass_filter(ts, cutoff, fs, order)
        norm_data[key] = ts

        if np.max(ts) > max_fluor:
            max_fluor = np.max(ts)

        ax.plot(time_vect,ts,color=cmap(col_indx/10))
        # ax.set_xlim(0,20000)
        col_indx+=1
    xl = ax.get_xlim()
    max_fluor = max_fluor/4
    # ax.plot([xl[0],xl[1]],[max_fluor,max_fluor],color='black',linewidth=2)
    ax.set_xlim(0,30000)
    row_indx+=1

# %%

# 2,3; 4,5; 6,7; 8,9; 10,11 are all technical replicates

col_pairs = [(2,3), (4,5), (6,7), (8,9), (10,11)]
data_avg = {}

fig,ax_list = plt.subplots(nrows=4,ncols=2,figsize=(8,10))
ax_list_t = ax_list.reshape(-1)

condition_max = []

row_indx = 0
for row in row_list:
    sample_indx = 0

    max_t = 0
    for cp in col_pairs:
        ax = ax_list_t[row_indx]
        key1 = row + str(cp[0])
        key2 = row + str(cp[1])
        
        ts_avg = (norm_data[key1] + norm_data[key2])/2

        if np.nanmax(ts_avg)>max_t:
            max_t = np.max(ts_avg)

        data_avg[row+str(sample_indx)] = ts_avg
        si = sample_times[sample_indx]*60
        ax.plot(time_vect-si,ts_avg,color=cmap(sample_indx/5))
        # ax.plot(time_vect-si,norm_data[key1],color=cmap(sample_indx/5))
        # ax.plot(time_vect-si,norm_data[key2],'--',color=cmap(sample_indx/5))
        # ax.set_xlim(0,10000)
        sample_indx+=1

    condition_max.append(max_t)
    row_indx+=1



# %% time to 250,000 RFU

# thresh = 200000
fig,ax_list = plt.subplots(nrows=4,ncols=2,figsize=(8,10))
ax_list_t = ax_list.reshape(-1)

# sample_times = [0,30,90,210,390]

time_to_thresh_data = {}

row_indx = 0
for row in row_list[0:-1]:
    ax = ax_list_t[row_indx]
    time_to_thresh = []
    thresh = condition_max[row_indx]/2

    for si in range(5):

        key = row + str(si)
        ts = data_avg[key]
        time_indx = np.argwhere(ts>=thresh)[0][0]
        time_t = time_vect[time_indx] - sample_times[si]*60
        time_to_thresh.append(time_t)

    ax.plot(sample_times,time_to_thresh)
    time_to_thresh_data[row] = time_to_thresh
    row_indx+=1


# %%

ydata0 = time_to_thresh_data['A']
ydata0 = ydata0 - np.min(ydata0)

growth_rates = []

for key in time_to_thresh_data.keys():
    ydata = time_to_thresh_data[key]

    # ydata = ydata-np.min(ydata)

    # ydata[0] = ydata0[0]
    st = sample_times
    ydata = ydata - ydata[0]
    # if key == 'E':
    #     ydata = np.delete(ydata,2)
        
    #     st = np.delete(st,2)

    popt = est_time_kill(st,ydata,debug=True)

    r = popt[0]
    A = popt[1]

    # if A > 0:
    #     r = -r
        

    growth_rates.append(-r)

# %%

growth_rates = growth_rates/np.max(growth_rates)

fig,ax = plt.subplots()
dc_plot = [-4,-3,-2,-1,0,1,2]
ax.scatter(dc_plot,growth_rates)
# ax.set_xscale('log')

# %%
