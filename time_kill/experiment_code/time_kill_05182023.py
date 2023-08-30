#%%
import re
import os
from fears.utils import AutoRate
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
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

with open('calibration_05172023/spline.pkl','rb') as f:
    spline = pickle.load(f)

with open('calibration_05172023/rfu30_to_cell_count.pkl','rb') as f:
    rfu30_to_cell_count = pickle.load(f)

exp_folder = 'tk_05182023/'
plate_paths = os.listdir(exp_folder)
min_cell_count = 3000

plate_paths = [p for p in plate_paths if p.endswith('.xlsx')]

plate_paths.sort(key=natural_keys)

drug_conc_labels = ['nd','10$^{-1}$','10$^{0}$','10$^{1}$','10$^{2}$','10$^{3}$']
#%%

plates = []
data = []
for plate_path in plate_paths:
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
            cell_count_t.append(rfu30_to_cell_count(data_t[key],spline))
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
fig,ax_list = plt.subplots(nrows=2,sharex=True)
# ax_list = ax_list.flatten()
cmap = mpl.colormaps['viridis']

row_indx = 0
sample_time = np.arange(0,10,1)*30

for row in row_list[1:-1]:
    if row_indx < 2:
        ax = ax_list[0]
    else:
        ax = ax_list[1]
    y = np.array(cell_count[row])
    st = sample_time
    yerr = np.array(cell_count_err[row])

    # ax.plot(fluor_data[row],label=row,color=cmap(row_indx/6))
    
    ax.errorbar(st,y,yerr=yerr,label=drug_conc_labels[row_indx],color=cmap(row_indx/5))
    row_indx += 1

    ax.set_yscale('log')

ax_list[0].set_ylabel('Cell count',fontsize=12)
ax_list[1].set_ylabel('Cell count',fontsize=12)
ax_list[1].set_xlabel('Sample time (min)',fontsize=12)

for ax in ax_list:
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

fig.legend(fontsize=10,frameon=False,bbox_to_anchor=(1.05,0.9),title='xMIC')
fig.suptitle('Raw cell count over time',fontsize=14)
# %%
fig,ax = plt.subplots()
# ax_list = ax_list.flatten()
cmap = mpl.colormaps['viridis']

row_indx = 0
sample_time = np.arange(0,10,1)*30

for row in row_list[1:-1]:

    y = np.array(cell_count[row])/cell_count[row][0]
    st = sample_time
    yerr = np.array(cell_count_err[row])/cell_count[row][0]
    # ax.plot(fluor_data[row],label=row,color=cmap(row_indx/6))
    
    ax.errorbar(st,y,yerr=yerr,label=drug_conc_labels[row_indx],color=cmap(row_indx/5))
    row_indx += 1

    ax.set_yscale('log')

ax.set_ylabel('Fold change in cells',fontsize=14)
ax.set_xlabel('Sample time (min)',fontsize=14)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

ax.legend(fontsize=10,frameon=False,title='xMIC')
# %%

# %%

cmap = mpl.colormaps['viridis']

def growth_diffeq(N,t,K,Kss,alpha,cc):
    # cc = 8.55 # carrying capacity
    # if N <= 0:
    #     dydt = 0
    # else:   
    dydt = (K-Kss*(1-np.exp(-alpha*t)))*N*(1-N/cc)

    return dydt

def growth_sol(t,y0,K,Kss,alpha,cc):
    y = odeint(growth_diffeq,y0,t,args=(K,Kss,alpha,cc))
    return y[:,0]

# def growth_sol_log(t,y0,K,Kss,alpha,cc):
#     y = odeint(growth_diffeq,y0,t,args=(K,Kss,alpha,cc))
#     return np.log(y[:,0])

cell_count_t = np.log10(cell_count['B']/cell_count['B'][0]) + 1
cell_count_err_t = np.log10(cell_count_err['B'])

p0 = [0,0.1,0,0,2]

st = sample_time
popt,pcov = curve_fit(growth_sol,st,cell_count_t,p0=p0)

# fig,ax_list = plt.subplots(nrows=2,sharex=True)
fig,ax = plt.subplots()
# ax = ax_list[0]
# ax.errorbar(sample_time[1:],cell_count['B'][1:],cell_count_err['B'][1:],label='data')
# ax.errorbar(st,cell_count_t,yerr=cell_count_err_t,color=cmap(0),fmt='o')
ax.scatter(st,cell_count_t,label='data',color=cmap(0))
xfit = np.linspace(0,300,100)
# popt = np.array([0.1, 0.03186479, 0.        , 0.        , 1.63826138])
yfit = growth_sol(xfit,*popt)
# yfit = (10**yfit)*cell_count['B'][1]
ax.plot(xfit,yfit,label='fit',color=cmap(0))

K = popt[1]
cc = popt[4]

st = sample_time
row_indx = 1

Kss = [0]
alpha = [0]
Kss_err = [0]

for row in row_list[2:-1]:
    # if row_indx < 2:
    #     ax = ax_list[0]
    # else:
    #     ax = ax_list[1]

    y = np.log10(np.array(cell_count[row])/cell_count[row][0]) + 1

    if row == 'F':
        y = y[0:4]    
        st = sample_time[0:4]
    elif row == 'G':
        y = y[0:3]
        st = sample_time[0:3]
    else:
        st = sample_time
   
    yerr = np.array(cell_count_err[row])
    
    p0 = [0,0.01,0]
    bounds = [[0,0,0],[5,0.1,np.inf]]
    popt,pcov = curve_fit(lambda st, y0, Kss, alpha: growth_sol(st,y0,K,Kss,alpha,cc),
                          st,y,p0=p0,maxfev=10000,bounds=bounds)
    row_indx += 1
    yfit = growth_sol(xfit,popt[0],K,popt[1],popt[2],cc)
    ax.plot(xfit,yfit,color=cmap(row_indx/6))
    ax.scatter(st,y,color=cmap(row_indx/6))
    # ax.set_yscale('log')
    Kss.append(popt[1])
    alpha.append(popt[2])
    Kss_err.append(np.sqrt(np.diag(pcov))[1])
# ax.set_ylim(0,5)
# %%
net_growth_rate = np.array(K-Kss)
# estimate growth rate for row G

cell_count_log = np.log10(cell_count['G'])
dy = np.diff(cell_count_log)
dt = np.diff(sample_time)
growth_rate = np.min(dy/dt)
net_growth_rate[-1] = growth_rate
# estimate growth rate for row F
cell_count_log = np.log10(cell_count['F'])
dy = np.diff(cell_count_log)
dt = np.diff(sample_time)
growth_rate = np.min(dy/dt)
net_growth_rate[-2] = growth_rate


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


fig,ax = plt.subplots()

drug_conc = [0,10**-1,10**0,10**1,10**2,10**3]
ax.scatter(drug_conc,net_growth_rate)
xt = ax.get_xticks()
# ax.set_xticks(drug_conc)
# ax.set_xticklabels(drug_conc_labels)
ax.set_xscale('log')

drug_conc = [0,10**-1,10**0,10**1,10**2,10**3]

max_gr = np.max(net_growth_rate)
min_gr = np.min(net_growth_rate)
mic = 1

p0 = [max_gr,min_gr,mic,.6]
bounds = [[max_gr - 0.1*max_gr, min_gr + 0.1*min_gr, mic - 0.1*mic,0],
          [max_gr + 0.1*max_gr, min_gr - 0.1*min_gr, mic + 0.1*mic,np.inf]]
popt,pcov = curve_fit(pharmacodynamic_curve,drug_conc,net_growth_rate,p0=p0,bounds=bounds)

xfit = np.logspace(-2,3,100)
yfit = pharmacodynamic_curve(xfit,*popt)

ax.plot(xfit,yfit)


# %%
