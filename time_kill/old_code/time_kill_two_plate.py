#%%
from fears.utils import AutoRate
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.signal import butter, lfilter
import scipy
from scipy.integrate import odeint
import pickle

# data_file = 'EK_AB_20230324_090141.xlsx'
data_file = 'EK_AB_20230330_101642.xlsx'

row_list = ['B','C','D','E','F','G']
col_list = [2,3,4,5,6,7,8,9,10,11]
col_list = [str(col) for col in col_list]

# sample_times = [0,30,90,210,390]
# # 11:12:44
# st_data_cleaning = [32*60-13,93*60+7,213*60+-18,397*60-26]

# 10:01:12
st_data_cleaning = [24*60+12,55*60+21,89*60+27,203*60-13]
sample_times = [0,26,56,89,233]

dc = [0,10**-1,10**0,10,10**2,10**3]
mic = 0.1

dc = [c*mic for c in dc]

filehandler = open('normalizer.p', 'rb') 
norm = pickle.load(filehandler)

#%% Helper functions
def get_piecewise_slope(xdata,ydata):
    slopes = []
    for i in range(len(xdata)-1):
        slopes.append((ydata[i]-ydata[i+1])/(xdata[i]-xdata[i+1]))
    return np.array(slopes)

def est_time_kill(xdata,ydata,debug=False):
    # alpha,A,l
    
    indx = np.argwhere(np.abs(ydata)==np.max(np.abs(ydata)))[0][0]
    A_est = ydata[indx]

    alpha_est = 0
    p0 = [alpha_est,A_est,0]
    # bounds = [[-np.inf,A_est-0.2*(np.abs(A_est)),-np.inf],
    #     [np.inf,A_est+0.2*(np.abs(A_est)),np.inf]]

    popt,pcov = scipy.optimize.curve_fit(time_kill_model,xdata,ydata,p0=p0,maxfev=10000)

    if debug:

        y_opt = time_kill_model(xdata,popt[0],popt[1],popt[2])
        
        fig,ax = plt.subplots()
        ax.plot(xdata,y_opt,label='est')
        ax.scatter(xdata,ydata,marker="*",label='data')
        ax.set_title(popt[0])
        ax.legend()
        # ax.set_ylim(-1000,50000)

    return popt

def net_growth(t,N0,K,Kss,alpha):

    g = []

    for tt in t:
        g.append(N0 + (1/2.303)*((K-Kss)*tt + (Kss/alpha)*(1-np.exp(-alpha*tt))))

    return g

def est_net_growth(xdata,ydata,K=None,N0=None,debug=False,p0=None,prev_alpha=None,alpha_diff=0.1,
                   min_alpha_diff=None):
    
    if K is None: # want to estimate K
        if p0 is None:
            p0 = [ydata[0],0.001,0,0.001] # N0, K, Kss, alpha
        bounds = [[ydata[0]-0.1*ydata[0],0,0,0.001],
                  [ydata[0]+0.1*ydata[0],np.inf,np.inf,1]]

        popt,pcov = scipy.optimize.curve_fit(lambda t,N0,K,Kss,alpha: net_growth(t,N0,K,Kss,alpha),
                                    xdata,ydata,p0 = p0,bounds=bounds,maxfev=10000)
    else:
        if p0 is None:
            if prev_alpha is not None:
                p0 = [0,prev_alpha + prev_alpha*alpha_diff] # Kss, alpha
            else:
                p0 = [0,0.01]
        if prev_alpha is not None:
            if min_alpha_diff is None:
                min_alpha_diff = alpha_diff/10
            bounds = [[0,prev_alpha + prev_alpha*min_alpha_diff],[np.inf,prev_alpha*1000]]
        else:
            bounds = [[0,10**-6],[10,1]]
        popt,pcov = scipy.optimize.curve_fit(lambda t,Kss,alpha: net_growth(t,N0,K,Kss,alpha),
                                    xdata,ydata,p0 = p0,maxfev=10000,bounds=bounds)
    
    if debug:

        if K is None:
            Kss = popt[2]
            y_opt = net_growth(xdata,popt[0],popt[1],popt[2],popt[3])
        else:
            Kss = popt[0]
            y_opt = net_growth(xdata,N0,K,popt[0],popt[1])
        
        fig,ax = plt.subplots()
        ax.plot(xdata,y_opt,label='est')
        ax.scatter(xdata,ydata,marker="*",label='data')
        ax.set_title(round(Kss,3))
        ax.legend()

    return popt
    

def time_kill_model(t,alpha,A,l):
    res = []
    for tt in t:
        res.append(A*(1-np.exp(-alpha*(tt-l)/A)))
    return res

# def time_kill_model(t,kss,alpha,y0):
#     # Growth rate constant of control arbitrarily set to zero
#     res = []
#     for tt in t:
#         res.append(y0 - kss*tt + (kss/alpha)*(1-np.exp(-1*alpha*tt)))
#     return res

def butter_lowpass(cutoff, fs, order=5):
    return butter(order, cutoff, fs=fs, btype='low', analog=False)

def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y

order = 1
fs = 0.01      # sample rate, Hz
cutoff = 0.0001  # desired cutoff frequency of the filter, Hz

def est_pharm_curve(xdata,ydata,debug=False):


    # gmax,gmin,mic,k

    g_max_est = ydata[0]
    gmin_est = np.min(ydata)

    mic_est = np.argwhere(np.array(ydata)<=0)[0][0]
    mic_est = xdata[mic_est]

    p0 = [g_max_est,gmin_est,mic_est,1]
    bounds = [[g_max_est-0.1*g_max_est,gmin_est+0.1*gmin_est,mic_est/2,0.1],
              [g_max_est+0.1*g_max_est,gmin_est-0.1*gmin_est,mic_est*2,10]]

    # p0 = [1,g_drugless_est,-1,0]
    # bounds = ([-np.inf,g_drugless_est-g_drugless_est*0.05,-2,0],
    #           [np.inf,g_drugless_est+g_drugless_est*0.05,-0.1,np.inf])

    popt,pcov = scipy.optimize.curve_fit(pharmacodynamic_curve,xdata,ydata,p0=p0,
                                         bounds=bounds)

    if debug:

        # xmin = np.log10(xdata[1])
        # xmax = np.log10(xdata[-1])
        dc_fit = np.logspace(-3,3)
        dc_fit[0] = 0
        y_opt = pharmacodynamic_curve(dc_fit,popt[0],popt[1],popt[2],popt[3])
        fig,ax = plt.subplots()
        ax.plot(dc_fit,y_opt)
        ax.plot(xdata,ydata,'o')

        ax.set_xscale('log')

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
        if c_t == 0:
            g.append(gmax)
        else:
            g.append(gmax - (((gmax-gmin)*(c_t/mic)**k)/((c_t/mic)**k-(gmin/gmax))))
    
    return g

# def modified_hill(c,)

#%%

col_pairs = [(2,3), (4,5), (6,7), (8,9), (10,11)]
data_avg = {}
data = {}

p = AutoRate.Plate(data_file)
time_vect = np.array(p.data['Time [s]'])

data_cleaning_indx = []
for tt in st_data_cleaning:
    indx = np.argwhere(time_vect>=tt)[0][0]
    p.data = p.data.drop(p.data.index[indx-1:indx+2])
    data_cleaning_indx.append(indx)
    time_vect = np.array(p.data['Time [s]'])

p.data = p.data.drop(p.data.index[-1]) # remove last datapoint because of incomplete data
time_vect = np.array(p.data['Time [s]'])

fig,ax_list = plt.subplots(nrows=3,ncols=2,figsize=(8,10))
ax_list_t = ax_list.reshape(-1)
cmap = mpl.colormaps['viridis']

row_indx = 0
for row in row_list:
    ax = ax_list_t[row_indx]
    sample_indx = 0

    for cp in col_pairs:
        key1 = row + str(cp[0])
        key2 = row + str(cp[1])

        ts1 = np.array(p.data[key1])
        ts1[ts1=='OVER'] = np.nan

        ts2 = np.array(p.data[key2])
        ts2[ts2=='OVER'] = np.nan
        
        ts_avg = np.nanmean((ts1,ts2),axis=0)

        data_avg[row+str(sample_indx)] = ts_avg
        data[row+str(sample_indx)] = (ts1,ts2)
        si = sample_times[sample_indx]*60

        # ts1 = butter_lowpass_filter(ts1, cutoff, fs, order)
        # ts2 = butter_lowpass_filter(ts2, cutoff, fs, order)
        # ax.plot(time_vect-si,ts_avg,color=cmap(sample_indx/5))
        ax.plot(time_vect-si,ts_avg,color=cmap(sample_indx/5))
        # ax.plot(time_vect-si,ts2,'--',color=cmap(sample_indx/5))
        ax.set_ylim(0,100000)
        ax.set_xlim(0,3600)
        sample_indx+=1
        ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0))

    row_indx+=1

#%% fluorescence at t = 3600

time_indx = np.argwhere(time_vect>=3600)[0][0]
# fig,ax_list = plt.subplots(nrows=3,ncols=2,figsize=(8,10))
# ax_list_t = ax_list.reshape(-1)

fig2,ax_all = plt.subplots()


# def rfu_to_cell_count(rfu):
#     logRFU = np.log10(rfu)
#     logN = 3.52*np.log10(21*(logRFU-4.02))

#     return 10**logN

def rfu_to_cell_count(rfu):

    res = []

    for r in rfu:
        dilution_factor = 8.66 + (-1.812*10**-5)*r
        dilution_factor=2**(-dilution_factor)
        res.append(dilution_factor*90000)
    return res

cell_count = []

cal_norm_factor = 1

fluor_3600 = []
samples = ['0','1','2','3','4']
row_indx = 0
for row in row_list:
    fluor_t = []
    sample_indx = 0
    for s in samples:
        key = row + s
        ts = data_avg[key]

        si = sample_times[sample_indx]*60
        time_indx = np.argwhere(time_vect-si>=3600)[0][0]
        # if dc[row_indx] > 0:
        #     f = ts[time_indx]/norm.get_ctx_norm_factor(np.log10(dc[row_indx]),3600)
        # else:
        f = ts[time_indx]
        fluor_t.append(f*cal_norm_factor)
        sample_indx+=1
    # ax = ax_list_t[row_indx]
    cc_t = rfu_to_cell_count(fluor_t)

    cell_count.append(cc_t)

    # ax.plot(sample_times,np.log10(cc_t))
    ax_all.plot(sample_times,cc_t,color=cmap(row_indx/5),label=round(dc[row_indx],2))
    # ax.set_ylim(4.05,5.65)
    row_indx+=1

ax_all.set_yscale('log')
ax_all.set_ylabel('Cell count (per uL)',fontsize=12)
ax_all.set_xlabel('Time (min)',fontsize=12)
ax_all.legend(frameon=False)
# ax_all.set_ylim(100,10000)

fig.savefig('time_kill_zoomed.png',bbox_inches='tight')
#%%

sample_times_hr = np.array(sample_times)/60

def growth_diffeq(N,t,K,Kss,alpha,cc):
    # cc = 8.55 # carrying capacity
    if N <= 0:
        dydt = 0
    else:   
        dydt = (K-Kss*(1-np.exp(-alpha*t)))*N*(1-N/cc)

    return dydt

def growth_sol(t,y0,K,Kss,alpha,cc):
    y = odeint(growth_diffeq,y0,t,args=(K,Kss,alpha,cc))
    return y[:,0]

cc0 = cell_count[0]

p0 = [10**3,1,0.01,10**6]
bounds = [[0.9*10**3,0,0,0],[10**4,10,1,10**7]]

popt,pcov = scipy.optimize.curve_fit(lambda t,y0,K,alpha,cc: growth_sol(t,y0,K,0,alpha,cc),
                                     sample_times_hr,cc0,p0=p0,bounds=bounds,maxfev=10000)

yt = growth_sol(sample_times_hr,popt[0],popt[1],0,popt[2],popt[3])

fig,ax_list = plt.subplots(nrows=6,figsize=(3,8))
ax = ax_list[0]
ax.plot(sample_times_hr,yt)
ax.scatter(sample_times_hr,cc0)
ax.set_yscale('log')
ax.set_ylim(500,5*10**5)
ax.set_xlim(-0.1,4.1)

ax.set_title('K = ' + str(round(popt[1],2)))

K = popt[1]
cc = popt[3]

p0 = [10**3,0,1]
bounds = [[0.5*10**3,0,0.5],[10**4,10,10]]

alpha_est = []
kss_est = []
kss_est.append(0)
alpha_est.append(0)

indx = 1
for c_count in cell_count[1:]:
    ax = ax_list[indx]
    st = sample_times_hr[c_count>500]
    c_count = c_count[c_count>500]
    popt,pcov = scipy.optimize.curve_fit(lambda t,y0,Kss,alpha: growth_sol(t,y0,K,Kss,alpha,cc),
                                     st,c_count,p0=p0,bounds=bounds,maxfev=10000)
    
    yt = growth_sol(st,popt[0],K,popt[1],popt[2],cc)
    ax.plot(st,yt)
    ax.scatter(st,c_count)
    ax.set_yscale('log')
    ax.set_ylim(500,5*10**5)
    ax.set_xlim(-0.1,4.1)
    alpha_est.append(popt[2])
    kss_est.append(popt[1])

    ax.set_title(r'$\alpha$ = ' + str(round(popt[2],2)) + ', $R\'_{ss}$ = ' + str(round(popt[1],2)))

    indx+=1
fig.tight_layout()
#%%

# # %% time to 200,000 

# thresh = 190000
# res = {}

# samples = ['0','1','2','3','4']

# fig,ax_list = plt.subplots(nrows=3,ncols=2,figsize=(8,10))
# ax_list_t = ax_list.reshape(-1)

# row_indx = 0
# for row in row_list:

#     time_to_thresh = []

#     sample_indx = 0

#     for s in samples:
#         key = row + s

#         ts = data_avg[key]
        
#         indx = np.argwhere(ts>=thresh)[0][0]

#         time_t = time_vect[indx] - sample_times[sample_indx]*60

#         time_to_thresh.append(time_t)
#         sample_indx +=1

#     ax = ax_list_t[row_indx]
#     ax.plot(sample_times,time_to_thresh)
#     res[str(row_indx)] = time_to_thresh

#     row_indx+=1

# # %%

# growth_rates = []

# fig,ax = plt.subplots()
# # ax_list_t = ax_list.reshape(-1)

# y0 = res['0'][0]
# cell_count = []

# indx = 0

# for key in res.keys():

#     ydata = np.array(res[key])

#     ydata[0] = y0

#     # convert to cell count
#     ydata = (-1.7*10**-4)*ydata - 2.96*10**-1
#     # ydata = 10**ydata

#     cell_count.append(ydata + 8)

#     st = sample_times

#     # if key == '4':
#     #     ydata = ydata[0:-1]
#     #     st = st[0:-1]

#     # if key == '3':
#     #     ydata = ydata[0:-1]
#     #     st = st[0:-1]

#     ax.plot(st,ydata,color=cmap(indx/5),label=round(dc[indx],2))

#     ydata = ydata - ydata[0]

#     indx+=1

#     popt = est_time_kill(st,ydata,debug=False)

#     r = popt[0]
#     A = popt[1]


#     growth_rates.append(r)

# ax.set_xlabel('Time (min)',fontsize=14)
# ax.set_ylabel('Log proportion of carrying capacity',fontsize=14)
# ax.legend(frameon=False)

# # fig,ax = plt.subplots()
# # dc_plot = [-2,-1,0,1,2,3]
# # ax.scatter(dc_plot,growth_rates)
# # ax.set_xscale('log')
# # # %%

# # popt = est_pharm_curve(dc,growth_rates,debug=True)

# # fig,ax = plt.subplots()

# # dc_fit = np.linspace(-3,3,num=100)
# # dc_fit = 10**dc_fit
# # dc_fit[0] = 0
# # g_fit = pharmacodynamic_curve(dc_fit,popt[0],popt[1],popt[2],1)
# # dc_fit_plot = np.linspace(-3,3,num=100)

# # ax.plot(dc_fit_plot,g_fit)
# # %%

# cell_count_trunc = []
# sample_times_trunc = []

# for i in range(len(cell_count)):
    
#     # if i in [0,1,3,4]:
#     #     cell_count_trunc.append(cell_count[i][0:-1])
#     #     sample_times_trunc.append(sample_times[0:-1])
#     # # elif i == 5:
#     # #     cell_count_trunc.append(cell_count[i][0:-2])
#     # #     sample_times_trunc.append(sample_times[0:-2])
#     # else:
#     cell_count_trunc.append(cell_count[i])
#     sample_times_trunc.append(sample_times)

# # killing_rate = []
# # alpha_est = []
# # alpha_fit = []

# # # estimate K
# # st = sample_times[0:len(cell_count_trunc[0])]
# # popt = est_net_growth(st,cell_count_trunc[0],debug=True)
# # K = popt[1]
# # # K = 0.008515791380952379
# # N0 = popt[0]

# # killing_rate.append(popt[2])
# # alpha_est.append(popt[-1])

# # alpha_fit.append(popt[-1])

# # indx = 0
# # for cc in cell_count_trunc[1:]:
# #     st = sample_times[0:len(cc)]
# #     popt = est_net_growth(st,cc,K=K,N0=N0,debug=True,prev_alpha=alpha_fit[-1],alpha_diff=1)
# #     # popt = est_net_growth(st,cc,K=K,N0=N0,debug=True)
# #     killing_rate.append(popt[0])
# #     if indx == 0:
# #         alpha_fit.append(popt[-1]/100)
# #     else:
# #         alpha_fit.append(popt[-1])
# #     alpha_est.append(popt[-1])
# #     indx += 1


# # %% Try ODEint

# sample_times_hr = np.array([t/60 for t in sample_times])

# def growth_diffeq(logN,t,K,Kss,alpha,cc):
#     # cc = 8.55 # carrying capacity
#     if logN <= 0:
#         dydt = 0
#     else:   
#         dydt = (K-Kss*(1-np.exp(-alpha*t)))*logN*(1-logN/cc)

#     return dydt

# def growth_sol(t,y0,K,Kss,alpha,cc):
#     y = odeint(growth_diffeq,y0,t,args=(K,Kss,alpha,cc))
#     return y[:,0]


# p0 = [6.5,0.1,0,0.2,7.5]
# bounds = [[0,0,0,10**-1,6],[10,10,10,0.3,8]]
# popt,pcov = scipy.optimize.curve_fit(growth_sol,sample_times_hr,cell_count[0],maxfev=10000,p0=p0,bounds=bounds)
# K = popt[1]
# N0 = popt[0]
# carry_cap = popt[4]

# fig,ax_list = plt.subplots(nrows=6,figsize=(4,11))
# ax = ax_list[0]
# yest = growth_sol(sample_times_hr,popt[0],popt[1],popt[2],popt[3],popt[4])
# ax.plot(sample_times_hr,yest)
# ax.scatter(sample_times_hr,cell_count[0])

# alpha_est = []
# killing_rate = []

# alpha_est.append(popt[3])
# killing_rate.append(popt[2])

# prev_alpha = 0.1
# indx = 1
# for c_count in cell_count[1:]:
# # for c_count in [cell_count_trunc[2]]:

#     ax = ax_list[indx]
#     st = sample_times_hr[c_count>0]
#     c_count = c_count[c_count>0]

#     # get Kss upper bound
#     slopes = get_piecewise_slope(st,c_count)

#     if any(slopes<0):
#         max_Kss = -np.min(slopes)
#     else:
#         max_Kss = 10

#     p0 = [0,1.1]
#     bounds = [[0,0.1],[max_Kss+0.1*max_Kss,100]]
    
#     popt,pcov = scipy.optimize.curve_fit(lambda t,Kss,alpha: growth_sol(t,N0,K,Kss,alpha,carry_cap),
#                                          st,c_count,maxfev=10000,p0=p0,bounds=bounds)

#     yest = growth_sol(st,N0,K,popt[0],popt[1],carry_cap)
#     ax.plot(st,yest)
#     ax.scatter(st,c_count)

#     alpha_est.append(popt[1])
#     killing_rate.append(popt[0])
#     prev_alpha = popt[1]
#     indx+=1
#     # ax.set_ylim(4,8)
# # %%
# killing_rate_t = np.array(killing_rate)
# dc_t = np.array(dc)
# # killing_rate_t = killing_rate_t[[0,1,2,4,5]]
# # dc_t = dc_t[[0,1,2,4,5]]

# net_growth = np.array(K-killing_rate_t)

# popt = est_pharm_curve(dc_t,net_growth,debug=True)

# fig,ax = plt.subplots()

# dc_fit = np.linspace(-3,3,num=100)
# dc_fit = 10**dc_fit
# dc_fit[0] = 0
# g_fit = pharmacodynamic_curve(dc_fit,popt[0],popt[1],popt[2],popt[3])
# dc_fit_plot = np.linspace(-3,3,num=100)

# ax.scatter(dc_t,net_growth)
# ax.plot(dc_fit,g_fit)
# ax.set_xscale('log')
# # %% Same analysis but don't average time traces
# thresh = 190000
# res = {}

# samples = ['0','1','2','3','4']

# fig,ax_list = plt.subplots(nrows=3,ncols=2,figsize=(8,10))
# ax_list_t = ax_list.reshape(-1)

# row_indx = 0
# for row in row_list:

#     time_to_thresh = []

#     sample_indx = 0

#     for s in samples:
#         key = row + s

#         ts1 = data[key][0]
#         ts2 = data[key][1]
        
#         indx1 = np.argwhere(ts1>=thresh)[0][0]

#         time_t1 = time_vect[indx1] - sample_times[sample_indx]*60

#         indx2 = np.argwhere(ts2>=thresh)[0][0]

#         time_t2 = time_vect[indx2] - sample_times[sample_indx]*60

#         time_to_thresh.append(np.mean((time_t1,time_t2)))
#         sample_indx +=1

#     ax = ax_list_t[row_indx]
#     ax.plot(sample_times,time_to_thresh)
#     res[str(row_indx)] = time_to_thresh

#     row_indx+=1

# #%%
# fig,ax = plt.subplots()
# # ax_list_t = ax_list.reshape(-1)

# y0 = res['0'][0]
# cell_count = []
# # growth_rates = []

# indx = 0

# for key in res.keys():

#     ydata = np.array(res[key])

#     ydata[0] = y0

#     # convert to cell count
#     ydata = (-1.7*10**-4)*ydata - 2.96*10**-1
#     # ydata = 10**ydata

#     cell_count.append(ydata + 8)

#     st = sample_times

#     # if key == '4':
#     #     ydata = ydata[0:-1]
#     #     st = st[0:-1]

#     # if key == '3':
#     #     ydata = ydata[0:-1]
#     #     st = st[0:-1]

#     ax.plot(st,ydata,color=cmap(indx/5),label=round(dc[indx],2))

#     ydata = ydata - ydata[0]

#     indx+=1

#     # popt = est_time_kill(st,ydata,debug=False)

#     # r = popt[0]
#     # A = popt[1]


#     # growth_rates.append(r)

# ax.set_xlabel('Time (min)',fontsize=14)
# ax.set_ylabel('Log proportion of carrying capacity',fontsize=14)
# ax.legend(frameon=False)


# # %%
# p0 = [6.5,0.1,0,0.2,7.5]
# bounds = [[0,0,0,10**-1,6],[10,10,10,0.3,8]]
# popt,pcov = scipy.optimize.curve_fit(growth_sol,sample_times_hr,cell_count[0],maxfev=10000,p0=p0,bounds=bounds)
# K = popt[1]
# N0 = popt[0]
# carry_cap = popt[4]

# fig,ax_list = plt.subplots(nrows=6,figsize=(4,11))
# ax = ax_list[0]
# yest = growth_sol(sample_times_hr,popt[0],popt[1],popt[2],popt[3],popt[4])
# ax.plot(sample_times_hr,yest)
# ax.scatter(sample_times_hr,cell_count[0])

# alpha_est = []
# killing_rate = []

# alpha_est.append(popt[3])
# killing_rate.append(popt[2])

# prev_alpha = 0.1
# indx = 1
# for c_count in cell_count[1:]:
# # for c_count in [cell_count_trunc[2]]:

#     ax = ax_list[indx]
#     st = sample_times_hr[c_count>0]
#     c_count = c_count[c_count>0]

#     # get Kss upper bound
#     slopes = get_piecewise_slope(st,c_count)

#     if any(slopes<0):
#         max_Kss = -np.min(slopes)
#     else:
#         max_Kss = 10

#     p0 = [0,1.1]
#     bounds = [[0,0.1],[max_Kss+0.1*max_Kss,100]]
    
#     popt,pcov = scipy.optimize.curve_fit(lambda t,Kss,alpha: growth_sol(t,N0,K,Kss,alpha,carry_cap),
#                                          st,c_count,maxfev=10000,p0=p0,bounds=bounds)

#     yest = growth_sol(st,N0,K,popt[0],popt[1],carry_cap)
#     ax.plot(st,yest)
#     ax.scatter(st,c_count)

#     alpha_est.append(popt[1])
#     killing_rate.append(popt[0])
#     prev_alpha = popt[1]
#     indx+=1
#     # ax.set_ylim(4,8)
# # %%
# killing_rate_t = np.array(killing_rate)
# dc_t = np.array(dc)
# killing_rate_t = killing_rate_t[[0,1,2,3,5]]
# dc_t = dc_t[[0,1,2,3,5]]

# net_growth = np.array(K-killing_rate_t)
# net_growth = net_growth/(cell_count[0][0])

# popt = est_pharm_curve(dc_t,net_growth,debug=True)

# fig,ax = plt.subplots()

# dc_fit = np.linspace(-3,3,num=100)
# dc_fit = 10**dc_fit
# dc_fit[0] = 0
# g_fit = pharmacodynamic_curve(dc_fit,popt[0],popt[1],popt[2],popt[3])
# dc_fit_plot = np.linspace(-3,3,num=100)

# ax.scatter(dc_t,net_growth)
# ax.scatter(dc,K-killing_rate)
# ax.plot(dc_fit,g_fit)
# ax.set_xscale('log')
# # %% try with raw cell count instead of log

# def growth_diffeq(N,t,K,Kss,alpha,cc):
#     # cc = 8.55 # carrying capacity
#     if N <= 0:
#         dydt = 0
#     else:   
#         dydt = (K-Kss*(1-np.exp(-alpha*t)))*N*(1-N/cc)

#     return dydt

# def growth_sol(t,y0,K,Kss,alpha,cc):
#     y = odeint(growth_diffeq,y0,t,args=(K,Kss,alpha,cc))
#     return y[:,0]

# cell_count_raw = []
# for cc in cell_count:
#     cell_count_raw.append(10**cc)

# c_count = cell_count_raw[0]

# # norm_const = np.max(c_count)

# # c_count = c_count/norm_const

# p0 = [10**6,0.1,0,0.2,10**8]
# bounds = [[0,0,0,10**-1,0.8],[10**7,10,10,0.3,10**10]]
# popt,pcov = scipy.optimize.curve_fit(growth_sol,sample_times_hr,c_count,maxfev=10000,p0=p0,bounds=bounds)

# K = popt[1]
# N0 = popt[0]
# carry_cap = popt[4]

# fig,ax_list = plt.subplots(nrows=6,figsize=(4,11))
# ax = ax_list[0]

# yest = growth_sol(sample_times_hr,popt[0],popt[1],popt[2],popt[3],popt[4])
# ax.plot(sample_times_hr,yest)
# ax.scatter(sample_times_hr,c_count)

# alpha_est = []
# killing_rate = []

# alpha_est.append(popt[3])
# killing_rate.append(popt[2])

# prev_alpha = 0.1
# indx = 1
# for c_count in cell_count_raw[1:]:
# # for c_count in [cell_count_trunc[2]]:
    
#     ax = ax_list[indx]
#     st = sample_times_hr
#     # c_count = c_count[c_count>1]

#     # get Kss upper bound
#     slopes = get_piecewise_slope(st,c_count)

#     if any(slopes<0):
#         max_Kss = -np.min(slopes) + K
#     else:
#         max_Kss = 10

#     p0 = [0,1.1]
#     # bounds = [[0,0.5],[100,100]]
    
#     popt,pcov = scipy.optimize.curve_fit(lambda t,Kss,alpha: growth_sol(t,N0,K,Kss,alpha,carry_cap),
#                                          st,c_count,maxfev=10000,p0=p0)

#     yest = growth_sol(st,N0,K,popt[0],popt[1],carry_cap)
#     ax.plot(st,yest)
#     ax.scatter(st,c_count)

#     alpha_est.append(popt[1])
#     killing_rate.append(popt[0])
#     prev_alpha = popt[1]
#     indx+=1




# %%
