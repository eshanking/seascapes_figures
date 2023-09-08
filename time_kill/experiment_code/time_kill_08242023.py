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

with open('../rfu_to_dilution.pkl','rb') as f:
    rfu30_to_dilution = pickle.load(f)

exp_folder = '../experiment_data/tk_08242023/'
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
fluor_err = {'B':[],'C':[],'D':[],'E':[],'F':[],'G':[]}
cell_count = {'B':[],'C':[],'D':[],'E':[],'F':[],'G':[]}
cell_count_err = {'B':[],'C':[],'D':[],'E':[],'F':[],'G':[]}
row_list = ['A','B','C','D','E','F','G','H']
dc = {'A':0,'B':'nd','C':-1,'D':0,'E':1,'F':2,'G':3,'H':0}

for sample_num in range(10):
    sample_col = np.mod(3*sample_num,9) + 2
    # plate_num = int(np.floor(sample_num/3))
    data_t = data[sample_num]

    col_list = [sample_col,sample_col+1,sample_col+2]
    col_list = [str(c) for c in col_list]
    # print((plate_num,col_list))
    row_indx = 0
    for row in row_list[1:-1]:
        # dc_t = dc[row]
        fluor_t = []
        cell_count_t = []
        # sf = scale_factor[row_indx]
        for col in col_list:
            # print((sample_num,row,col))
            key = row + col
            fluor_t.append(data_t[key])
            # cell_count_t.append(rfu30_to_cell_count(data_t[key],spline))
            cell_count_t.append(rfu30_to_dilution(data_t[key]))
        row_indx += 1
        fluor_avg = np.mean(fluor_t)
        fluor_err_t = np.std(fluor_t)/np.sqrt(len(fluor_t))
        cell_count_avg = np.mean(cell_count_t)
        cell_count_err_t = np.std(cell_count_t)/np.sqrt(len(cell_count_t))
        fluor_data[row].append(fluor_avg)
        fluor_err[row].append(fluor_err_t)
        cell_count[row].append(cell_count_avg)
        cell_count_err[row].append(cell_count_err_t)

#%% plot raw fluorescence data

fig,ax_list = plt.subplots(nrows=2)
ax_list = ax_list.flatten()
cmap = mpl.colormaps['viridis']

row_indx = 0
sample_time = np.arange(0,10,1)*30

for row in row_list[1:-1]:
    if row_indx < 2:
        ax = ax_list[0]
    else:

        ax = ax_list[1]
    y = np.array(fluor_data[row])
    st = sample_time
    yerr = np.array(fluor_err[row])
    # ax.plot(fluor_data[row],label=row,color=cmap(row_indx/6))
    
    ax.errorbar(st,y,yerr=yerr,color=cmap(row_indx/5),label=drug_conc_labels[row_indx])
    row_indx += 1

    ax.set_yscale('log')

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

for row in row_list[1:-1]:
    # if row_indx < 2:
    #     ax = ax_list[0]
    # else:
    #     ax = ax_list[1]
    y = np.array(cell_count[row])
    y = np.log2(y)
    st = sample_time
    yerr = np.array(cell_count_err[row])
    # yerr = np.log2(yerr)

    # ax.plot(fluor_data[row],label=row,color=cmap(row_indx/6))
    y = y-y[0]
    ax.errorbar(st,y,yerr=yerr,label=drug_conc_labels[row_indx],color=cmap(row_indx/5))
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

def growth_diffeq(N,t,K,Kss,alpha,cc):
    # cc = 8.55 # carrying capacity
    # if N <= 0:
    #     dydt = 0
    # else:   
    dydt = (K-Kss*(1-np.exp(-alpha*t)))*N*(1-N/cc)

    return dydt

def growth_sol(t,K,Kss,alpha,cc):
    y = odeint(growth_diffeq,1,t,args=(K,Kss,alpha,cc))
    return y[:,0]

# estimate K

y = np.log10(cell_count['B'])
y = y-y[0] + 1

p0 = [0.01,0,0,4]

x = sample_time
popt,pcov = curve_fit(growth_sol,x,y,p0=p0)

# plot fit

fig,ax = plt.subplots()

ax.plot(x,y,'o',color='k',label='data')

xfit = np.linspace(0,300,100)
yfit = growth_sol(xfit,*popt)
ax.plot(xfit,yfit,'-',label='fit',color=cmap(0))

K = popt[0]
cc = popt[3]
row_indx = 1

alpha_list = [0]
Kss_list = [0]

prev_alpha = 0

for row in row_list[2:-1]:
    x = sample_time
    y = np.log10(cell_count[row])
    y = y-y[0] + 1

    if row == 'G':
        y = np.delete(y,1)
        x = np.delete(x,1)
        prev_alpha = 0.1
    if row == 'C':
        y = np.delete(y,np.s_[8:])
        x = np.delete(x,np.s_[8:])

    ax.plot(x,y,'o',color=cmap(row_indx/5),label=row)

    p0 = [0.1,0.1]
    bounds = [[0,prev_alpha],[1,1]]

    popt,pcov = curve_fit(lambda x, Kss, alpha: growth_sol(x,K,Kss,alpha,cc),
                        x,y,p0=p0,maxfev=10000,bounds=bounds)
    
    Kss = popt[0]
    alpha = popt[1]
    yfit = growth_sol(xfit,K,Kss,alpha,cc)
    ax.plot(xfit,yfit,'-',color=cmap(row_indx/5))

    alpha_list.append(alpha)
    Kss_list.append(Kss)
    prev_alpha = alpha

    row_indx += 1

# %%
net_rate = K - Kss_list

fig,ax = plt.subplots()

drug_conc = [0,0.01,0.1,1,10,100]

ax.plot(drug_conc,net_rate,'o',color='k')
ax.set_xscale('symlog',linthresh=0.01)
# %%
