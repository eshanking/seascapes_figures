#%%
import re
import os
from fears.utils import AutoRate
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.integrate import odeint
from scipy.optimize import curve_fit

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

def inv_logistic(y,y0,k,x50,l):
    return x50/((y0/(y-l) - 1)**k)

def RFU_to_cell_count(RFU):
    """Generates a cell count estimate from RFU data at 1 hour timepoint.

    RFU data was calculated with a constant gain of 40.

    Args:
        RFU (float): Measured RFU value from 1 hour timepoint

    Returns:
        float: Estimated cell count (cells/uL).
    """
     
    RFU_norm_const = 5868

    if RFU > RFU_norm_const:
        print('Warning: RFU value ' + str(round(RFU,2)) + ' is above calibration range.')
    if RFU < 500:
        # print('Warning: RFU value is below calibration range.')
        print('Warning: RFU value ' + str(round(RFU,2)) + ' is below calibration range.')
    
    cell_count_norm_const = 905401
    popt = [0.94254949, 0.77916359, 0.09317131, 0.10172957]

    est = inv_logistic(RFU/RFU_norm_const,*popt)*cell_count_norm_const

    return est

exp_folder = 'tk_05112023/'
plate_paths = os.listdir(exp_folder)

plate_paths = [p for p in plate_paths if p.endswith('.xlsx')]

plate_paths.sort(key=natural_keys)

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

sf = 2.5

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
            cell_count_t.append(RFU_to_cell_count(data_t[key]/sf))
        row_indx += 1
        fluor_avg = np.mean(fluor_t)
        fluor_err_t = np.std(fluor_t)/np.sqrt(len(fluor_t))
        cell_count_avg = np.mean(cell_count_t)
        cell_count_err_t = np.std(cell_count_t)/np.sqrt(len(cell_count_t))
        fluor_data[row].append(fluor_avg)
        fluor_err[row].append(fluor_err_t)
        cell_count[row].append(cell_count_avg)
        cell_count_err[row].append(cell_count_err_t)

# %%
fig,ax = plt.subplots()
ax_list = ax_list.flatten()
cmap = mpl.colormaps['viridis']

row_indx = 0
sample_time = np.arange(0,10,1)*30

for row in row_list[1:-1]:
    if row_indx < 2:
        # ax = ax_list[0] 
        y = np.array(cell_count[row][1:])/cell_count[row][1]
        st = sample_time[1:]
        yerr = np.array(cell_count_err[row][1:])/cell_count[row][1]
    else:
        # ax = ax_list[1]
        y = np.array(cell_count[row])/cell_count[row][0]
        st = sample_time
        yerr = np.array(cell_count_err[row])/cell_count[row][0]
    # ax.plot(fluor_data[row],label=row,color=cmap(row_indx/6))
    
    ax.errorbar(st,y,yerr=yerr,label=row,color=cmap(row_indx/6))
    row_indx += 1

    ax.set_yscale('log')

ax.set_ylabel('Fold change',fontsize=14)
ax.set_xlabel('Sample time (min)',fontsize=14)
# %%
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
    # ax = ax_list[row_indx]
    # ax.plot(fluor_data[row],label=row,color=cmap(row_indx/6))
    ax.errorbar(sample_time,fluor_data[row],yerr=fluor_err[row],label=row,color=cmap(row_indx/6))
    row_indx += 1

    # ax.set_yscale('log')
# %%
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

cell_count_norm = np.log10(cell_count['B'][1:])
cell_count_norm_err = cell_count_err['B'][1:]/cell_count['B'][1]

p0 = [0,0.1,0,0,2]

st = sample_time[1:]-sample_time[1]
popt,pcov = curve_fit(growth_sol,st,cell_count_norm,p0=p0)

fig,ax = plt.subplots()
# ax.errorbar(sample_time[1:],cell_count['B'][1:],cell_count_err['B'][1:],label='data')
ax.errorbar(st,cell_count_norm,label='data')
xfit = np.linspace(0,300,100)
# popt = np.array([0.1, 0.03186479, 0.        , 0.        , 1.63826138])
yfit = growth_sol(xfit,*popt)
# yfit = (10**yfit)*cell_count['B'][1]
ax.plot(xfit,yfit,label='fit')
# %%
