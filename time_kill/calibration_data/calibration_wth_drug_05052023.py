#%%
from fears.utils import AutoRate
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import stats
import scipy.interpolate as interp
import scipy.optimize as sciopt

ab_plate_path = 'calibration_05052023/EK_AB_constant_gain_20230505_145012.xlsx'
od_plate_path = 'calibration_05052023/EK_single_OD600_20230505_102953.xlsx'
p_ab = AutoRate.Plate(ab_plate_path)
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

dilutions = []
for i in range(len(col_list)-2):
    dilutions.append(str(2**i) + 'x')

ax.errorbar(dilutions,od_avg,yerr=od_std,fmt='o')

def od_to_cells(od):
    res = [1350069.89751786,   -7125.64194304]
    return res[1] + res[0]*od

cell_count = od_to_cells(od_avg)
cell_count = cell_count - np.min(cell_count)

time_vect = np.array(p_ab.data['Time [s]'])

# %%

fig,ax_list = plt.subplots(nrows=3,ncols=2)

ax_list = ax_list.flatten()
cmap = mpl.colormaps['viridis']

col_indx = 0

ax_indx = 0
for row in row_list[1:-1]:
    col_indx = 0
    for col in col_list[2:-1]:
        ax = ax_list[ax_indx]
        key = row + col
        ax.plot(time_vect/60,np.array(p_ab.data[key]).astype('float64'),color=cmap(col_indx/8))
        col_indx += 1
    ax_indx += 1

# for ax in ax_list:
    # ax.set_xlim(0,60)
    # ax.set_ylim(0,400000)

fig.tight_layout()

# %%
# fig,ax_list = plt.subplots(nrows=4,figsize=(4,8))
cmap = mpl.colormaps['viridis']

dilution = []
for i in range(len(col_list)-2):
    dilution.append(1/(2**i))

hr1_time_indx = np.argwhere(time_vect>=3600)[0][0]

col_indx = 0

RFU_60 = np.zeros((6,9))

col_indx = 0
for col in col_list[2:-1]:
    ts_t = np.zeros(len(time_vect))
    row_indx = 0
    for row in row_list[1:-1]:
        key = row + col
        ts = np.array(p_ab.data[key]).astype('float64')
        # no_drug_fluor.append(ts[hr1_time_indx])
        RFU_60[row_indx][col_indx] = ts[hr1_time_indx]
        row_indx += 1

    col_indx += 1

RFU60_avg = np.mean(RFU_60,axis=0)
RFU60_err = np.std(RFU_60,axis=0)/np.sqrt(6)

fig,ax = plt.subplots()

for row in range(6):
    ax.plot(cell_count[1:],RFU_60[row,:],color=cmap(row/5))

ax.set_xscale('log')

# %%
fig,ax = plt.subplots()
# cell_count[-1] = 1
ax.errorbar(RFU60_avg,cell_count[1:],xerr=RFU60_err,fmt='o',color='k',capsize=5)

def expon(x,y0,alpha,l):
    return y0*np.exp(alpha*(x-l))

cell_count_norm = cell_count[1:]/np.max(cell_count[1:])
RFU60_avg_norm = RFU60_avg/np.max(RFU60_avg)

ax.set_yscale('log')


def inv_logistic(y,y0,k,x50,l):
    return x50/((y0/(y-l) - 1)**k)

p0 = [1.1,1,0.1,0.1]
popt_weighted,pcov = sciopt.curve_fit(inv_logistic,RFU60_avg_norm,cell_count_norm,p0=p0,sigma=RFU60_err/np.max(RFU60_avg))
# popt,pcov = sciopt.curve_fit(expon,RFU60_avg,cell_count[1:],p0=[1,1,1])

xfit = np.linspace(0.01,1,100)
yfit = inv_logistic(xfit,*popt_weighted)

ax.plot(xfit*np.max(RFU60_avg),yfit*np.max(cell_count[1:]),label='weighted LS',linewidth=2)

popt,pcov = sciopt.curve_fit(inv_logistic,RFU60_avg_norm,cell_count_norm,p0=p0)
# popt,pcov = sciopt.curve_fit(expon,RFU60_avg,cell_count[1:],p0=[1,1,1])

xfit = np.linspace(0.01,1,100)
yfit = inv_logistic(xfit,*popt)

ax.plot(xfit*np.max(RFU60_avg),yfit*np.max(cell_count[1:]),'--',label='unweighted LS')

ax.legend(frameon=False)

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

ax.set_xlabel('RFU$_{60}$',fontsize=14)
ax.set_ylabel('Cell Count (cells/$\mu$L)',fontsize=14)
ax.tick_params(labelsize=12)
# %%

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

    if RFU > RFU_norm_const or RFU < 500:
        print('Warning: RFU value is outside of calibration range.')
    
    cell_count_norm_const = 905401
    popt = [0.94254949, 0.77916359, 0.09317131, 0.10172957]

    est = inv_logistic(RFU/RFU_norm_const,*popt)*cell_count_norm_const

    return est

# %%
