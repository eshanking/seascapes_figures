#%%
import re
import os
from fears.utils import AutoRate
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy
from scipy.integrate import odeint
from scipy.optimize import curve_fit
import pickle
import pandas as pd
import calibration_08302023 as cal

rfu30_to_dilution = cal.run()

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

# import pickled objects

# with open('calibration_05172023/spline.pkl','rb') as f:
#     spline = pickle.load(f)

# with open('../rfu_to_dilution.pkl','rb') as f:
#     rfu30_to_dilution = pickle.load(f)

exp_folder = '../experiment_data/tk_08302023'
plate_paths = os.listdir(exp_folder)
min_cell_count = 3000

plate_paths = [p for p in plate_paths if p.endswith('.xlsx')]

plate_paths.sort(key=natural_keys)

drug_conc_labels = ['nd','10$^{-1}$','10$^{0}$','10$^{1}$','10$^{2}$','10$^{3}$']
#%%

plates = []
data = []
for plate_path in plate_paths:
    if plate_path[0] != '~':
        path_t = os.getcwd() + os.sep + exp_folder + os.sep + plate_path
        p = AutoRate.Plate(path_t,mode='single_measurement')
        plates.append(p)
        data.append(p.od_data_to_dict(p.data))

# %%
# sample_num = 0
# sample_col = np.mod(3*sample_num,9) + 2

fluor_data = {'B':[],'C':[],'D':[],'E':[],'F':[],'G':[]}
# fluor_err = {'B':[],'C':[],'D':[],'E':[],'F':[],'G':[]}
cell_count = {'B':[],'C':[],'D':[],'E':[],'F':[],'G':[]}
# cell_count_err = {'B':[],'C':[],'D':[],'E':[],'F':[],'G':[]}
row_list = ['A','B','C','D','E','F','G','H']
dc = {'A':0,'B':'nd','C':-1,'D':0,'E':1,'F':2,'G':3,'H':0}

for sample_num in range(10):
    # plate_num = int(np.floor(sample_num/3))
    data_t = data[sample_num]

    sample_col = str(sample_num + 2)
    # print((plate_num,col_list))
    row_indx = 0
    for row in row_list[1:-1]:

        fluor_t = []
        cell_count_t = []
        key = row + sample_col
        fluor = data_t[key]
        cell_count_t = rfu30_to_dilution(fluor)
        row_indx += 1
        fluor_data[row].append(fluor)
        cell_count[row].append(cell_count_t)

#%% Remove time point 3 because lid was left on

for row in row_list[1:-1]:
    fluor_data[row] = np.delete(fluor_data[row],2)
    cell_count[row] = np.delete(cell_count[row],2)


#%% plot raw fluorescence data

fig,ax_list = plt.subplots(nrows=2)
ax_list = ax_list.flatten()
cmap = mpl.colormaps['viridis']

row_indx = 0
sample_time = np.arange(0,10,1)*30
sample_time = np.delete(sample_time,2)

for row in row_list[1:-1]:
    if row_indx < 2:
        ax = ax_list[0]
    else:

        ax = ax_list[1]
    y = np.array(fluor_data[row])
    st = sample_time
    # ax.plot(fluor_data[row],label=row,color=cmap(row_indx/6))
    
    ax.plot(st,y,color=cmap(row_indx/5),label=drug_conc_labels[row_indx])
    row_indx += 1

    # ax.set_yscale('log')

ax_list[0].set_ylabel('Fluorescence (RFU)',fontsize=12)
ax_list[1].set_ylabel('Fluorescence (RFU)',fontsize=12)
ax.set_xlabel('Sample time (min)',fontsize=12)

# set fig title
fig.suptitle('Raw fluorescence over time',fontsize=14)

for ax in ax_list:
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

fig.legend(fontsize=10,frameon=False,bbox_to_anchor=(1.05,0.9),title='xMIC')

#%% Plot raw cell count data
# fig,ax_list = plt.subplots(nrows=2,sharex=True)
fig,ax = plt.subplots(figsize=(6,4))
# ax_list = ax_list.flatten()
cmap = mpl.colormaps['viridis']

row_indx = 0
sample_time = np.arange(0,10,1)*30
sample_time = np.delete(sample_time,2)

for row in row_list[1:-1]:
    # if row_indx < 2:
    #     ax = ax_list[0]
    # else:
    #     ax = ax_list[1]
    y = np.array(cell_count[row])
    y = np.log2(y)
    st = sample_time
    # yerr = np.log2(yerr)

    # ax.plot(fluor_data[row],label=row,color=cmap(row_indx/6))
    y = y-y[0]
    ax.plot(st,y,label=drug_conc_labels[row_indx],color=cmap(row_indx/5))
    row_indx += 1

    # ax.set_yscale('log')

# ax_list[0].set_ylabel('Log$_{2}$ cell count',fontsize=12)
# ax_list[1].set_ylabel('Cell count',fontsize=12)
# ax_list[1].set_xlabel('Sample time (min)',fontsize=12)

# for ax in ax_list:

ax.set_ylabel('Log$_{2}$ fold change',fontsize=12)

ax.set_xticks(sample_time)
ax.set_xticklabels(sample_time/60)

ax.set_xlabel('Sample time (hr)',fontsize=12)

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

fig.legend(fontsize=10,frameon=False,bbox_to_anchor=(1.05,0.9),title='xMIC')
fig.suptitle('Cell count over time',fontsize=14)

# %%
cmap = mpl.colormaps['viridis']

def growth_diffeq(N,t,K,Kss,alpha,cc):
    # cc = 8.55 # carrying capacity
    # if N <= 0:
    #     dydt = 0
    # else:   
    dydt = (K-Kss*(1-np.exp(-alpha*t)))*N*(1-N/cc)

    return dydt

def growth_sol(t,K,Kss,alpha,cc):
    y = odeint(growth_diffeq,2,t,args=(K,Kss,alpha,cc))
    return y[:,0]

# estimate K

y = np.log10(cell_count['B'])
y = y-y[0] + 2

p0 = [0.01,4]
alpha = 0
Kss = 0

x = sample_time
x = np.array(x,dtype="float64")
y = np.array(y,dtype="float64")

rate_err = []

popt,pcov = curve_fit(lambda x, K, cc: growth_sol(x,K,Kss,alpha,cc),x,y,p0=p0,
                      maxfev=10000)

K = popt[0]
cc = popt[1]

rate_err.append(np.sqrt(np.diag(pcov))[0])

# plot fit

fig,ax = plt.subplots()

ax.plot(x,y,'o',color='k',label='data')

xfit = np.linspace(0,np.max(x),100)
yfit = growth_sol(xfit,K,0,0,cc)
ax.plot(xfit,yfit,'-',label='fit',color=cmap(0))

# K = popt[0]
# cc = popt[3]
row_indx = 1

alpha_list = [0]
Kss_list = [0]

prev_alpha = 0

for row in row_list[2:-1]:
    x = sample_time
    y = np.log10(cell_count[row])
    y = y-y[0] + 2

    #     prev_alpha = 0.1
    # if row == 'C':
    #     y = np.delete(y,np.s_[8:])
    #     x = np.delete(x,np.s_[8:])

    ax.plot(x,y,'o',color=cmap(row_indx/5),label=row)

    p0 = [0.1,0.1]
    bounds = [[0,0],[1,1]]
    # bounds = [[0,prev_alpha],[1,1]]

    popt,pcov = curve_fit(lambda x, Kss, alpha: growth_sol(x,K,Kss,alpha,cc),
                        x,y,p0=p0,maxfev=10000,bounds=bounds)
    
    rate_err.append(np.sqrt(np.diag(pcov))[0])

    Kss = popt[0]
    alpha = popt[1]
    yfit = growth_sol(xfit,K,Kss,alpha,cc)
    ax.plot(xfit,yfit,'-',color=cmap(row_indx/5))

    alpha_list.append(alpha)
    Kss_list.append(Kss)
    prev_alpha = alpha

    row_indx += 1

net_rate = K - np.array(Kss_list)

rate_err = np.array(rate_err)
#%%

# def growth_diffeq(N,t,K,Kss,alpha,cc):
#     # cc = 8.55 # carrying capacity
#     # if N <= 0:
#     #     dydt = 0
#     # else:   
#     dydt = (K-Kss*(1-np.exp(-alpha*t)))*N*(1-N/cc)

#     return dydt

# def growth_sol(t,K,Kss,alpha,cc):
#     y = odeint(growth_diffeq,1,t,args=(K,Kss,alpha,cc))
#     return y[:,0]

# # estimate K

# y = np.log10(cell_count['B'])
# y = y-y[0] + 1

# p0 = [0.01,0,0,4]

# x = sample_time
# popt,pcov = curve_fit(growth_sol,x,y,p0=p0)

# # plot fit

# fig,ax = plt.subplots()

# ax.plot(x,y,'o',color='k',label='data')

# xfit = np.linspace(0,300,100)
# yfit = growth_sol(xfit,*popt)
# ax.plot(xfit,yfit,'-',label='fit',color=cmap(0))

# K = popt[0]
# cc = popt[3]
# row_indx = 1

# alpha_list = [0]
# Kss_list = [0]

# prev_alpha = 0

# for row in row_list[2:-1]:
#     x = sample_time
#     y = np.log10(cell_count[row])
#     y = y-y[0] + 1

#     if row == 'G':
#         y = np.delete(y,1)
#         x = np.delete(x,1)
    
#     x = np.delete(x,y<0)
#     y = np.delete(y,y<0)
#     #     prev_alpha = 0.1
#     # if row == 'C':
#     #     y = np.delete(y,np.s_[8:])
#     #     x = np.delete(x,np.s_[8:])

#     ax.plot(x,y,'o',color=cmap(row_indx/5),label=row)

#     p0 = [0.1,0.1]
#     bounds = [[0,0],[1,1]]
#     # bounds = [[0,prev_alpha],[1,1]]

#     popt,pcov = curve_fit(lambda x, Kss, alpha: growth_sol(x,K,Kss,alpha,cc),
#                         x,y,p0=p0,maxfev=10000,bounds=bounds)
    
#     Kss = popt[0]
#     alpha = popt[1]
#     yfit = growth_sol(xfit,K,Kss,alpha,cc)
#     ax.plot(xfit,yfit,'-',color=cmap(row_indx/5))

#     alpha_list.append(alpha)
#     Kss_list.append(Kss)
#     prev_alpha = alpha

#     row_indx += 1


# net_rate = K - Kss_list

#%%

fig,ax = plt.subplots()

drug_conc = [0,0.01,0.1,1,10,100]

ax.plot(drug_conc,net_rate,'o',color='k')
ax.set_xscale('symlog',linthresh=0.01)
# %% manually choose linear range

sub_mic_range = [0,1,2]
above_mic_range = [1,2,3,4,5,6,7,8]

fig,ax = plt.subplots()

rate = []

rate_err = []

row_indx = 0
for row in row_list[1:-1]:
    if row == 'B' or row == 'C':
        range_t = sub_mic_range
        
    else:
        range_t = above_mic_range
    
    x = sample_time[range_t]
    y = np.log10(cell_count[row]*10**6)[range_t]
    
    ax.plot(x,y,'o',color=cmap(row_indx/5),label=row)

    res = scipy.stats.linregress(x,y)
    rate.append(res.slope)
    rate_err.append(res.stderr)

    # plot fit
    xfit = np.linspace(0,300,100)
    yfit = res.slope*xfit + res.intercept
    ax.plot(xfit,yfit,'-',color=cmap(row_indx/5))

    row_indx += 1
rate = np.array(rate)
rate_err = np.array(rate_err)
# %% fit pharmacodynamic model

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

gmax = rate[0]

gmin = np.min(rate)

# fit mic and k

p0 = [0.1,1]
popt,pcov = curve_fit(lambda x, mic, k: pharmacodynamic_curve(x,gmax,gmin,mic,k),
                        drug_conc[1:],rate[1:],p0=p0,maxfev=10000)

mic = popt[0]
k = popt[1]

xfit = np.logspace(-3,2,100)
yfit = np.array(pharmacodynamic_curve(xfit,gmax,gmin,mic,k))

fig,ax = plt.subplots()

drug_conc = [0,0.01,0.1,1,10,100]

# ax.plot(drug_conc,rate*60,'o',color='k')

ax.errorbar(drug_conc,rate*60,yerr=rate_err*60,fmt='o',color='k',
            label='data')

ax.plot(xfit,yfit*60,'-',color='red')

ax.set_xscale('symlog',linthresh=0.01)

df = pd.DataFrame({'drug_conc':drug_conc,'rate':rate*60,'mic':mic,'k':k,
                   'err':rate_err*60})

df.to_csv('results_08312023.csv',index=False)

# %%
