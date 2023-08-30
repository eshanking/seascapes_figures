#%%
from fears.utils import AutoRate
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import stats
import scipy.interpolate as interp
import scipy.optimize as sciopt

ab_plate_path = 'calibration_05162023/EK_single_AB_constant_gain_20230516_114142.xlsx'
od_plate_path = 'calibration_05162023/EK_single_OD600_20230516_103116_with_lid.xlsx'
p_ab = AutoRate.Plate(ab_plate_path,mode='single_measurement')
p_od = AutoRate.Plate(od_plate_path,mode='single_measurement')
od_data = p_od.od_data_to_dict(p_od.data)

#%% estimate od background

row_list = ['A','B','C','D','E','F','G','H']
bg_cols = [1,12]
bg_cols = [str(c) for c in bg_cols]

bg_est = 0
indx = 0
for row in row_list:
    for col in bg_cols:
        key = row + col
        bg_est += od_data[key]
        indx+=1

bg_est = bg_est/indx

col_list = np.arange(12) + 1
col_list = [str(c) for c in col_list]

for row in row_list:
    for col in col_list:
        key = row+col
        od_data[key] = od_data[key] - bg_est

#%% Estimate cell count from OD
row_indx = 0

od_data_t = np.zeros((len(row_list)-2,len(col_list)-2))
for row in row_list[1:-1]:
    col_indx = 0
    for col in col_list[1:-1]:
        key = row + col
        od_data_t[row_indx,col_indx] = od_data[key]
        col_indx += 1
    row_indx += 1

od_avg = np.mean(od_data_t,axis=0)
od_std = np.std(od_data_t,axis=0)

fig,ax = plt.subplots()

dilutions_str = []
dilutions = []
for i in range(len(col_list)-2):
    dilutions_str.append(str(2**i) + 'x')
    dilutions.append(1/(2**i))

ax.errorbar(dilutions_str,od_avg,yerr=od_std,fmt='x',capsize=5,color='k')

def od_to_cells(od):
    res = 1350069.89751786
    return res*od

cell_count = od_to_cells(od_avg)

cell_count_est = [d*cell_count[0] for d in dilutions]

fig,ax = plt.subplots()

ax.plot(cell_count_est,cell_count,'o',color='k')
ax.set_xscale('log')
ax.set_yscale('log')

(ymin,ymax) = ax.get_ylim()
(xmin,xmax) = ax.get_xlim()

ax.set_ylim([np.min([xmin,ymin]),np.max([xmax,ymax])])  
ax.set_xlim([np.min([xmin,ymin]),np.max([xmax,ymax])]) 

ax.set_xlabel('Dilution estimation method',fontsize=14)
ax.set_ylabel('OD600 estimation method',fontsize=14)

# linear fit
res = stats.linregress(cell_count_est,cell_count)

xfit = np.linspace(np.min(cell_count_est),np.max(cell_count_est),100)
yfit = res.slope*xfit + res.intercept
ax.plot(xfit,yfit,'--',color='k',label='linear fit')


cell_count = cell_count_est
# annotate r^2



# K = 9*10**5 # cells/uL

# cell_count = K*np.array(dilutions)

# %%
# fig,ax_list = plt.subplots(nrows=4,figsize=(4,8))
fluor_data = p_ab.od_data_to_dict(p_ab.data)

cmap = mpl.colormaps['viridis']

fluor_data_t = np.zeros((len(row_list)-2,len(col_list)-2))
fluor_data_std = np.zeros((len(row_list)-2,len(col_list)-2))
row_indx = 0
for row in row_list[1:-1]:
    col_indx = 0
    for col in col_list[1:-1]:
        key = row + col
        fluor_data_t[row_indx,col_indx] = fluor_data[key]
        col_indx += 1
    row_indx += 1

fluor_avg = np.mean(fluor_data_t,axis=0)
fluor_err = np.std(fluor_data_t,axis=0)/np.sqrt(len(row_list)-2)

fig,ax = plt.subplots()

ax.errorbar(fluor_avg[1:],cell_count[1:],xerr=fluor_err[1:],fmt='o',capsize=5,color='k')

# ax.set_xscale('log')
ax.set_yscale('log')
#%%
fig,ax = plt.subplots()

ax.errorbar(fluor_avg[1:],cell_count[1:],xerr=fluor_err[1:],fmt='o',capsize=5,color='k')

cell_count_norm = cell_count[1:]/np.max(cell_count[1:])
fluor_avg_norm = fluor_avg[1:]/np.max(fluor_avg[1:])

def inv_logistic(y,y0,k,x50,l):
    return x50/((y0/(y-l) - 1)**k)

def power_law_func(X, a, b, c, d):
    # X = X - l
    return a * np.power(X, b) + c * np.power(X, d)

p0 = [1.1,1,0.1,0.1]
popt,pcov = sciopt.curve_fit(power_law_func,fluor_avg_norm,cell_count_norm,
                                      p0=p0,sigma=fluor_err[1:]/np.max(fluor_avg[1:]))
# # popt,pcov = sciopt.curve_fit(expon,RFU60_avg,cell_count[1:],p0=[1,1,1])
# popt[-1] = 
xfit = np.linspace(np.min(cell_count_norm),1,100)
yfit = power_law_func(xfit,*popt)

ax.plot(xfit*np.max(fluor_avg[1:]),yfit*np.max(cell_count[1:]),label='weighted LS',linewidth=2)

ax.set_yscale('log')

# %%