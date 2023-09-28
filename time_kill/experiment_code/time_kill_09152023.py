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

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

# with open('../rfu_to_dilution.pkl','rb') as f:
#     rfu30_to_dilution = pickle.load(f)
#     f.close()

rfu30_to_dilution = cal.run()

plate1_path = '../experiment_data/tk_09152023/plate_1'
plate2_path = '../experiment_data/tk_09152023/plate_2'

drug_conc = [0,10,100,1000] # x MIC
drug_conc = np.array(drug_conc)*0.1

plate_layout = {0:(1,['E','F','G']),1:(1,['B','C','D']),
                10:(2,['E','F','G']),100:(2,['B','C','D'])} # first index is plate number, second index is rows

sample_time = [0,15,35,70,113,170,235,295,355]

def get_timeseries(plate_path):
    plate_files = os.listdir(plate_path)
    plate_files = [p for p in plate_files if p.endswith('.xlsx')]
    plate_files.sort(key=natural_keys)
    plates = []
    data = []
    for plate_file in plate_files:
        if plate_file[0] != '~':
            path_t = os.getcwd() + os.sep + plate_path + os.sep + plate_file
            p = AutoRate.Plate(path_t,mode='single_measurement')
            plates.append(p)
            data.append(p.od_data_to_dict(p.data))
    return plates,data

plates1,data1 = get_timeseries(plate1_path)
plates2,data2 = get_timeseries(plate2_path)

data = (data1,data2)
#%%

# plate 1 row E,F,G: no drug
# plate 1 row B,C,D: 10x MIC
# plate 2 row E,F,G: 100x MIC
# plate 2 row B,C,D: 1000x MIC

timeseries_dict = {}
mean_dict = {}
err_dict = {}

cell_count_dict = {}
cell_count_err_dict = {}

for key in plate_layout:
    plate_num = plate_layout[key][0]
    rows = plate_layout[key][1]
    data_t = np.zeros((len(rows),len(sample_time)))
    row_indx = 0
    for row in rows:
        col_indx = 0
        for col_num in range(len(sample_time)):
            col_num_str = str(col_num + 2) # account for plate indexing
            key_t = row + col_num_str
            data_t[row_indx,col_indx] = data[plate_num-1][col_indx][key_t]
            col_indx += 1
        row_indx += 1
    timeseries_dict[key] = data_t
    mean_dict[key] = np.mean(data_t,axis=0)
    err_dict[key] = np.std(data_t,axis=0)/np.sqrt(len(rows))

    cell_count_t = rfu30_to_dilution(data_t)
    cell_count_dict[key] = np.mean(cell_count_t,axis=0)
    cell_count_err_dict[key] = np.std(cell_count_t,axis=0)/np.sqrt(len(rows))


# %%

# fig,ax = plt.subplots()

# for key in mean_dict:
#     y = mean_dict[key]
#     y = y[1:]
#     y = y-y[0]
#     x = sample_time[1:]
#     yerr = err_dict[key][1:]
#     ax.errorbar(x,y,yerr=yerr,label=str(key) + 'x MIC')

# ax.set_xlabel('Time (min)',fontsize=14)

# ax.set_title('Fluorescence Data')
# ax.set_ylabel('Fluorescence (RFU)',fontsize=14)
# ax.legend(frameon=False)

# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)

# ax.ticklabel_format(style='sci',scilimits=(0,0),axis='y')

# %%

fig,ax = plt.subplots()

for key in mean_dict:
    y = cell_count_dict[key]
    y = y[1:]
    y = y-y[0]
    x = sample_time[1:]
    yerr = cell_count_err_dict[key][1:]
    ax.errorbar(x,y,yerr=yerr,label=str(key) + 'x MIC')

ax.set_xlabel('Time (min)',fontsize=14)

ax.set_title('Cell Count Data')
ax.set_ylabel('Fold Change in Cell Count',fontsize=14)
ax.legend(frameon=False)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

ax.ticklabel_format(style='sci',scilimits=(0,0),axis='y')


# %%
cmap = mpl.colormaps['viridis']

def growth_diffeq(N,t,K,Kss,alpha,cc):

    dydt = (K-Kss*(1-np.exp(-alpha*t)))*N*(1-N/cc)

    return dydt

def growth_sol(t,y0,K,Kss,alpha,cc):
    y = odeint(growth_diffeq,y0,t,args=(K,Kss,alpha,cc))
    return y[:,0]

# estimate K

y = np.log10(cell_count_dict[0])
y = y - y[0] + 2

p0 = [y[0],0.01,4]
alpha = 0
Kss = 0

x = sample_time

rate_err = []

popt,pcov = curve_fit(lambda x, y0, K, cc: growth_sol(x,y0,K,Kss,alpha,cc),x,y,p0=p0,
                      maxfev=10000)

y0 = popt[0]
K = popt[1]
cc = popt[2]

rate_err.append(np.sqrt(np.diag(pcov))[0])

# plot fit

fig,ax_list = plt.subplots(nrows=4,figsize=(4,10),sharex=True,sharey=False)

ax = ax_list[0]
# ax.plot(x,y,'o',color='k',label='data')
init_count = np.log10(cell_count_dict[0]*10**6)[0]
y = y + init_count
ax.errorbar(x,y,yerr=cell_count_err_dict[0],fmt='o',color='k')

xfit = np.linspace(0,np.max(x),100)
yfit = growth_sol(xfit,y0,K,0,0,cc)
yfit = yfit + init_count

ax.plot(xfit,yfit,'-',label='fit',color=cmap(0))


row_indx = 1

alpha_list = [0]
Kss_list = [0]

prev_alpha = 0

for dc in drug_conc[1:]:
    ax = ax_list[row_indx]
    x = sample_time
    y = np.log10(cell_count_dict[dc])
    y = y - y[0] + 2
    init_count = np.log10(cell_count_dict[dc]*10**6)[0]

    ax.plot(x,y + init_count,'o',color=cmap(row_indx/5),label=row)

    p0 = [2,0.1,0.1]
    bounds = [[1,0,0],[3,1,1]]
    # bounds = [[0,prev_alpha],[1,1]]

    popt,pcov = curve_fit(lambda x, y0, Kss, alpha: growth_sol(x,y0, K,Kss,alpha,cc),
                        x,y,p0=p0,maxfev=10000,bounds=bounds)
    
    rate_err.append(np.sqrt(np.diag(pcov))[0])

    y0 = popt[0]
    Kss = popt[1]
    alpha = popt[2]
    yfit = growth_sol(xfit,y0,K,Kss,alpha,cc)

    yfit = yfit + init_count

    ax.plot(xfit,yfit,'-',color=cmap(row_indx/5))

    alpha_list.append(alpha)
    Kss_list.append(Kss)
    prev_alpha = alpha

    row_indx += 1

net_rate = K - np.array(Kss_list)

rate_err = np.array(rate_err)

ax_list[-1].set_xlabel('Time (min)',fontsize=12)

row_indx = 0

for ax in ax_list:
    ax.set_ylabel('log$_{10}$ Cell Count',fontsize=12)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='both',labelsize=12)

    ax.annotate(str(drug_conc[row_indx]) + ' ug/mL',xy=(0.7,0.5),xycoords='axes fraction',fontsize=12)
    row_indx += 1


# %% Load previous results

df = pd.read_csv('results_08312023.csv')

prev_dc = df['drug_conc'].to_numpy()
prev_rate = df['rate'].to_numpy()
prev_err = df['err'].to_numpy()

prev_rate = np.array(prev_rate)/60

#%% fit pharmacodynamic curve

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

rate = np.array(net_rate)

gmax = rate[0]

gmin = np.min(rate)

# fit k

# mic = 0.1

p0 = [1,0.1]

dc_t = np.concatenate((drug_conc,prev_dc))
rate_t = np.concatenate((rate,prev_rate))

popt,pcov = curve_fit(lambda x, k, mic: pharmacodynamic_curve(x,gmax,gmin,mic,k),
                        dc_t,rate_t,p0=p0,maxfev=10000)

# mic = popt[0]
k = popt[0]
mic = popt[1]

xfit = np.logspace(-3,2,100)
xfit = np.concatenate(([0],xfit))
yfit = np.array(pharmacodynamic_curve(xfit,gmax,gmin,mic,k))

fig,ax = plt.subplots(figsize=(4,3))

ax.plot(xfit,yfit*60,'-',linewidth=3)

# ax.errorbar(drug_conc,rate*60,yerr=rate_err*60,color='k',
#             label='09152023',capsize=3,fmt="D")
# # ax.plot(prev_dc,prev_rate*60,'o',color='red',label='08312023')
# ax.errorbar(prev_dc,prev_rate*60,yerr=prev_err,fmt='o',color='red',
#             label='08312023',capsize=3)
ax.scatter(drug_conc,rate*60,color='k',
            label='09152023',marker="D")

ax.scatter(prev_dc,prev_rate*60,color='red',
            label='08312023')

ax.set_xscale('symlog',linthresh=0.001)

ax.legend(frameon=False)

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

ax.set_xlabel('Drug Concentration (ug/mL)',fontsize=14)
ax.set_ylabel('Net Growth Rate (h$^{-1}$)',fontsize=14)

df = pd.DataFrame({'g_drugless':gmax*60,'gmin':gmin*60,'mic':mic,'k':k},index=[0])

df.to_csv('../results/pharm_params_09152023.csv',index=False)

# ax.set_xscale('symlog',linthresh=0.01)

# %%
