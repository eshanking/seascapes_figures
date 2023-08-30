#%%
from fears.utils import AutoRate
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy
import scipy.interpolate as interp

ab_plate_path = "calibration_04132023/EK_AB_20230414_111621.xlsx"
od_plate_path = "calibration_04132023/EK_single_OD600_20230413_120140.xlsx"
p_ab = AutoRate.Plate(ab_plate_path)
p_od = AutoRate.Plate(od_plate_path, mode='single_measurement')
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

#%% estimate cell count per column

def od_to_cells(od):
    res = [297761.03865714, 324941.3815491, 17609.09483884]
    return res[2] + res[1]*od + res[0]*od**2

col_list = np.arange(10) + 2
row_list = ['B','C','D','E','F','G']
cell_count = []
od = []

for col in col_list:
    od_avg = 0
    for row in row_list:
        key = row + str(col)
        od_avg += od_data[key]
    od_avg = od_avg/6
    cell_count.append(od_to_cells(od_avg))
    od.append(od_avg)

cell_count_est = cell_count

cell_count_col_1 = cell_count[0]
cell_count = []

for i in range(10):
    cell_count.append(cell_count_col_1/(2**i))

# %%
time_vect = np.array(p_ab.data['Time [s]'])

fig,ax_list = plt.subplots(nrows=2)
cmap = mpl.colormaps['viridis']

col_indx = 0

for col in col_list[1:]:
    ax = ax_list[0] 
    ts_avg = np.zeros(len(time_vect))
    for row in ['B','C','D']:
        key = row + str(col)
        ts_avg += np.array(p_ab.data[key]).astype('float64')
    ts_avg = ts_avg/3
    ax.plot(time_vect,ts_avg,color=cmap(col_indx/8))
    # col_indx += 1
    ax = ax_list[1] 
    ts_avg = np.zeros(len(time_vect))
    for row in ['E','F','G']:
        key = row + str(col)
        ts_avg += np.array(p_ab.data[key]).astype('float64')
    ts_avg = ts_avg/3
    ax.plot(time_vect,ts_avg,color=cmap(col_indx/8))
    col_indx += 1

for ax in ax_list:
    ax.set_xlim(0,3600)
    ax.set_ylim(0,400000)
# %%
fig,ax_list = plt.subplots(nrows=4,figsize=(4,8))
cmap = mpl.colormaps['viridis']

dilution = []
for i in range(len(col_list)):
    dilution.append(str(2**i) + 'x')

col_indx = 0

dilutions_no_drug = np.zeros(9)
dilutions_drug = np.zeros((3, 9))
ind60 = np.where((time_vect/60)>60)[0][0]

for col in col_list[1:]:
    ax = ax_list[0] 
    ts_avg = np.zeros(len(time_vect))
    for row in ['B','C','D']:
        key = row + str(col)
        ts_avg += np.array(p_ab.data[key]).astype('float64')
    ts_avg = ts_avg/3
    dilutions_no_drug[col_indx] = ts_avg[ind60]
    ax.plot(time_vect/60,ts_avg,color=cmap(col_indx/8))

    row_indx = 1
    for row in ['E','F','G']:
        ax = ax_list[row_indx]
        key = row + str(col)
        ts = np.array(p_ab.data[key]).astype('float64')
        ax.plot(time_vect/60,ts,color=cmap(col_indx/8),label=dilution[col_indx])
        dilutions_drug[row_indx-1][col_indx] = ts[ind60]
        row_indx+=1
    col_indx+=1

ax_list[0].set_title('No drug')
ax_list[1].set_title('1x MIC')
ax_list[2].set_title('10x MIC')
ax_list[3].set_title('100x MIC')

for ax in ax_list:
    ax.set_xlim(0,60)
    ax.set_ylim(0,500000)
    ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
    ax.set_xlabel('Time (min)')
    ax.set_ylabel('RFU (au)')

ax_list[3].legend(ncol=4,frameon=False,loc=(-0.05,-1.35))

fig.tight_layout()


# %%

# Graphs for fluorescence values vs. dilution at 60 minutes
print('Graphs for fluorescence vs. dilution using fluorescence values at 60 minutes.')

fig2,ax_list2 = plt.subplots(nrows=4,figsize=(4,8))

dilutionarray = np.array([1, 1./2, 1/4., 1./8, 1./16, 1./32, 1./64, 1./128, 1./256])
ax_list2[0].plot(dilutionarray, dilutions_no_drug, marker='o', markersize=3)
ax_list2[1].plot(dilutionarray, dilutions_drug[0], marker='o', markersize=3)
ax_list2[2].plot(dilutionarray, dilutions_drug[1], marker='o', markersize=3)
ax_list2[3].plot(dilutionarray, dilutions_drug[2], marker='o', markersize=3)

ax_list2[0].set_title('Fluorescence vs. Dilution (No drug)')
ax_list2[1].set_title('Fluorescence vs. Dilution (1x MIC)')
ax_list2[2].set_title('Fluorescence vs. Dilution (10x MIC)')
ax_list2[3].set_title('Fluorescence vs. Dilution (100x MIC)')

for ax in ax_list2:
    ax.set_ylim(0,500000)
    ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
    ax.set_xlabel('Dilution')
    ax.set_ylabel('RFU (au)')

fig2.tight_layout()


# %%

# Graphs for fluorescence values vs. dilution using area under the curve from 0-60 min
print('Graphs for fluorescence values vs. dilution using area under the curve from 0-60 min\n')
fig3,ax_list3 = plt.subplots(nrows=4,figsize=(4,8))

col_indx=0

curve_area_no_drug = np.zeros(9)
curve_area_drug = np.zeros((3, 9))

for col in col_list[1:]:
    ts_avg = np.zeros(len(time_vect))
    for row in ['B','C','D']:
        key = row + str(col)
        ts_avg += np.array(p_ab.data[key]).astype('float64')
    ts_avg = ts_avg/3
    ts_spline = interp.InterpolatedUnivariateSpline(time_vect[:-1:]/60, ts_avg[:-1:])
    xfine = np.linspace(0, 60, 1000)
    #ax_list3[0].plot(xfine, ts_spline(xfine, 0), color=cmap(col_indx/8))
    curve_area_no_drug[col_indx] = ts_spline.integral(0, 60)

    row_indx = 1
    for row in ['E','F','G']:
        ax = ax_list[row_indx]
        key = row + str(col)
        ts = np.array(p_ab.data[key]).astype('float64')
        ts_drug_spline = interp.InterpolatedUnivariateSpline(time_vect[:-1:]/60, ts[:-1:])
        #ax_list3[row_indx].plot(xfine,ts_drug_spline(xfine, 0),color=cmap(col_indx/8),label=dilution[col_indx])
        curve_area_drug[row_indx-1][col_indx] = ts_drug_spline.integral(0, 60)
        row_indx+=1

    col_indx+=1

print('Area under the curve from 0 to 60 minutes (no drug):', curve_area_no_drug)
print('Area under the curve from 0 to 60 minutes (with drug):', curve_area_drug)
print(dilutionarray)
ax_list3[0].plot(dilutionarray, curve_area_no_drug, marker='o', markersize=3)
ax_list3[1].plot(dilutionarray, curve_area_drug[0], marker='o', markersize=3)
ax_list3[2].plot(dilutionarray, curve_area_drug[1], marker='o', markersize=3)
ax_list3[3].plot(dilutionarray, curve_area_drug[2], marker='o', markersize=3)

ax_list3[0].set_title('Fluorescence vs. Dilution (No drug)')
ax_list3[1].set_title('Fluorescence vs. Dilution (1x MIC)')
ax_list3[2].set_title('Fluorescence vs. Dilution (10x MIC)')
ax_list3[3].set_title('Fluorescence vs. Dilution (100x MIC)')

for ax in ax_list3:
    #ax.set_ylim(0,500000)
    ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
    ax.set_xlabel('Dilution')
    ax.set_ylabel('RFU (au) \n(area under curve)')

fig3.tight_layout()

# %%
fig,ax_list = plt.subplots(nrows=4,figsize=(4,8))
cmap = mpl.colormaps['viridis']

dilution = []
for i in range(len(col_list)):
    dilution.append(1/(2**i))

hr1_time_indx = np.argwhere(time_vect>=3600)[0][0]

col_indx = 0

for col in col_list[1:]:
    ax = ax_list[0] 
    ts_avg = np.zeros(len(time_vect))
    for row in ['B','C','D']:
        key = row + str(col)
        ts_avg += np.array(p_ab.data[key]).astype('float64')
    ts_avg = ts_avg/3
    
    ax.plot(time_vect/60,ts_avg,color=cmap(col_indx/8))

    row_indx = 1
    for row in ['E','F','G']:
        ax = ax_list[row_indx]
        key = row + str(col)
        ts = np.array(p_ab.data[key]).astype('float64')
        ax.plot(time_vect/60,ts,color=cmap(col_indx/8),label=dilution[col_indx])
        row_indx+=1
    col_indx+=1
# %%
