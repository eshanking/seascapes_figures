#%%
from fears.utils import AutoRate
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import stats
import scipy.interpolate as interp
import scipy.optimize as sciopt
import pickle

# ab_plate_path = 'calibration_05172023/EK_single_AB_constant_gain_20230517_161938_30_min.xlsx'
ab_plate_path = 'calibration_05172023/EK_single_AB_constant_gain_20230517_165026_60_min.xlsx'
od_plate_path = 'calibration_05172023/EK_single_OD600_20230517_154009_no_lid.xlsx'
p_ab = AutoRate.Plate(ab_plate_path,mode='single_measurement')
p_od = AutoRate.Plate(od_plate_path,mode='single_measurement')
od_data = p_od.od_data_to_dict(p_od.data)

#%% Load pickled od to cell count function

with open('od_cell_count_no_lid_spline.pkl','rb') as f:
    od_to_cell_count = pickle.load(f)

#%% background subtract od data

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

od_t = np.zeros((len(row_list)-2,len(col_list)-2))
cell_count_t = np.zeros((len(row_list)-2,len(col_list)-2))
for row in row_list[1:-1]:
    col_indx = 0
    for col in col_list[1:-1]:
        key = row + col
        cell_count_t[row_indx,col_indx] = od_to_cell_count(od_data[key])
        od_t[row_indx,col_indx] = od_data[key]
        col_indx += 1
    row_indx += 1

cell_count_avg = np.mean(cell_count_t,axis=0)
cell_count_std = np.std(cell_count_t,axis=0)
cell_count_err = cell_count_std/np.sqrt(len(row_list)-2)

od_avg = np.mean(od_t,axis=0)
od_std = np.std(od_t,axis=0)
od_err = od_std/np.sqrt(len(row_list)-2)

fig,ax = plt.subplots()

dilutions_str = []
dilutions = []
for i in range(len(col_list)-2):
    dilutions_str.append(str(2**i) + 'x')
    dilutions.append(1/(2**i))

ax.errorbar(dilutions[1:],cell_count_avg[1:],yerr=cell_count_err[1:],fmt='o',capsize=5,color='k')
ax.set_xscale('log')
ax.set_yscale('log')

fig,ax = plt.subplots()

ax.errorbar(dilutions[1:],od_avg[1:],yerr=od_err[1:],fmt='o',capsize=5,color='k')
ax.set_xscale('log')
ax.set_yscale('log')
# annotate r^2



# K = 9*10**5 # cells/uL

# cell_count = K*np.array(dilutions)

# %%
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

ax.errorbar(fluor_avg[1:],cell_count_avg[1:],xerr=fluor_err[1:],yerr=cell_count_err[1:],fmt='o',capsize=5,color='k')

# ax.set_xscale('log')
ax.set_yscale('log')
#%%
fig,ax = plt.subplots()

ax.errorbar(fluor_avg[1:],cell_count_avg[1:],xerr=fluor_err[1:],yerr=cell_count_err[1:],fmt='o',capsize=5,color='k')

cell_count_norm = cell_count_avg[1:]/np.max(cell_count_avg[1:])
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
xfit = np.linspace(np.min(fluor_avg_norm),1,100)
yfit = power_law_func(xfit,*popt)

ax.plot(xfit*np.max(fluor_avg[1:]),yfit*np.max(cell_count_avg[1:]),label='power law',linewidth=2)

ax.set_yscale('log')

opt_err = np.sqrt(np.diag(pcov))
print('power law parameter variance: ',opt_err)

# estimate residuals

residuals = cell_count_norm - power_law_func(fluor_avg_norm,*popt)
print('\nresiduals: ',residuals)

# calculate RSS
RSS = np.sum(residuals**2)
print('\nRSS: ',RSS)

# try with inverse logistic

p0 = [1.1, 0.1, 0.5, 0] # y0,k,x50,l
bounds = ([0,0,0,0],[10,10,1,1])
popt,pcov = sciopt.curve_fit(inv_logistic,fluor_avg_norm,cell_count_norm,
                                      p0=p0,bounds=bounds,sigma=fluor_err[1:]/np.max(fluor_avg[1:]))

yfit = inv_logistic(xfit,*popt)
ax.plot(xfit*np.max(fluor_avg[1:]),yfit*np.max(cell_count_avg[1:]),label='inverse logistic',linewidth=2)

ax.legend(frameon=False)

opt_err = np.sqrt(np.diag(pcov))

print('\ninv logistic parameter variance: ',opt_err)

# estimate residuals
residuals = cell_count_norm - inv_logistic(fluor_avg_norm,*popt)
print('\nresiduals: ',residuals)

# calculate RSS
RSS = np.sum(residuals**2)
print('\nRSS: ',RSS)

# %%