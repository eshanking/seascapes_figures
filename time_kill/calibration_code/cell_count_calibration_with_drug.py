#%%
from fears.utils import AutoRate
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import stats
import scipy.interpolate as interp
import scipy.optimize as sciopt

ab_plate_path = 'time_kill/calibration_04132023/EK_AB_20230414_111621.xlsx'
od_plate_path = 'time_kill/calibration_04132023/EK_single_OD600_20230413_120140.xlsx'
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

#%% estimate ab background

time_vect = np.array(p_ab.data['Time [s]'])
bg = np.zeros(len(time_vect))

indx = 0
for row in row_list:
    key = row + '11'
    bg += np.array(p_ab.data[key]).astype('float64')
    indx += 1
bg = bg/indx

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
# fig,ax_list = plt.subplots(nrows=4,figsize=(4,8))
cmap = mpl.colormaps['viridis']

dilution = []
for i in range(len(col_list)):
    dilution.append(1/(2**i))

hr1_time_indx = np.argwhere(time_vect>=3600)[0][0]

col_indx = 0

no_drug_fluor = np.zeros((3,9))
no_drug_fluor_avg = []

col_indx = 0
for col in col_list[1:]:
    ts_avg = np.zeros(len(time_vect))
    row_indx = 0
    for row in ['B','C','D']:
        key = row + str(col)
        ts = np.array(p_ab.data[key]).astype('float64')
        ts_avg += ts
        # no_drug_fluor.append(ts[hr1_time_indx])
        no_drug_fluor[row_indx][col_indx] = ts[hr1_time_indx]
        row_indx += 1

    ts_avg = ts_avg/3
    no_drug_fluor_avg.append(ts_avg[hr1_time_indx])
    col_indx += 1

#1xMIC

MICx1_fluor = []

for col in col_list[1:]:

    key = 'E' + str(col)
    ts = np.array(p_ab.data[key]).astype('float64')

    MICx1_fluor.append(ts[hr1_time_indx])

#10xMIC

MICx10_fluor = []

for col in col_list[1:]:

    key = 'F' + str(col)
    ts = np.array(p_ab.data[key]).astype('float64')

    MICx10_fluor.append(ts[hr1_time_indx])

#100xMIC

MICx100_fluor = []

for col in col_list[1:]:

    key = 'G' + str(col)
    ts = np.array(p_ab.data[key]).astype('float64')

    MICx100_fluor.append(ts[hr1_time_indx])

fig,ax = plt.subplots()

dilution_log = np.arange(10)
# ax.plot(dilution_log[1:],no_drug_fluor,color=cmap(0),linewidth=2,label='no drug')
ax.plot(dilution_log[1:],MICx1_fluor,color=cmap(0.33),linewidth=2,label='1xMIC')
ax.plot(dilution_log[1:],MICx10_fluor,color=cmap(0.66),linewidth=2,label='10xMIC')
ax.plot(dilution_log[1:],MICx100_fluor,color=cmap(0.99),linewidth=2,label='100xMIC')

dilution_xlabel = ['$2^{' + str(d) + '}$' for d in dilution_log]
ax.set_xticks(np.arange(1,10))
ax.set_xticklabels(dilution_xlabel[1:])

ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
ax.tick_params(axis='both',labelsize=12)

ax.set_ylabel('RFU (a.u.)',fontsize=14)
ax.set_xlabel('Dilution',fontsize=14)

ax.legend(frameon=False,fontsize=12)
fig.savefig('cell_count_vs_dilution.png',bbox_inches='tight')
# %% linear fit

# ydata = np.mean((no_drug_fluor,MICx100_fluor,MICx10_fluor,MICx1_fluor),axis=0)

# res = stats.linregress(dilution_log[1:],ydata)

fig,ax = plt.subplots()

K = (432/50)*10**5 # carrying capacity cells per uL
cell_count = [K/(2**d) for d in dilution_log]
cell_count_log = np.log(cell_count)

rfu_data = np.concatenate((no_drug_fluor,[MICx1_fluor],[MICx10_fluor],[MICx100_fluor]))
rfu_data_mean = np.mean(rfu_data,axis=0)
rfu_err = np.std(rfu_data,axis=0)/np.sqrt(6)

# dilution_xlabel = ['$2^{' + str(d) + '}$' for d in dilution_log]

ax.errorbar(rfu_data_mean[1:],cell_count[2:],xerr=rfu_err[1:],fmt='o',color='black')

# ax.set_xticks(np.arange(1,10))
# ax.set_xticklabels(dilution_xlabel[1:])

ax.ticklabel_format(axis='x',style='sci',scilimits=(0,0))
ax.tick_params(axis='both',labelsize=12)

ax.set_xlabel('RFU$_{60}$ (a.u.)',fontsize=14)
ax.set_ylabel('Cell count (uL$^{-1}$)',fontsize=14)

# fit exponential

# def expon(x,y0,alpha,l):
#     return y0*np.exp(alpha*(x-l))

def inv_logistic(y,y0,k,x50):
    return x50/((y0/y - 1)**k)

rfu_norm = rfu_data_mean/np.max(rfu_data_mean)
cell_count_norm = cell_count/np.max(cell_count)

p0 = [1.1,1,0.01]
bounds = [[1.001,0,0],[10,10,0.06]]

popt,pcov = sciopt.curve_fit(inv_logistic,rfu_norm[1:],cell_count_norm[2:])

xfit = np.arange(0,1,0.01)
yfit = inv_logistic(xfit,popt[0],popt[1],popt[2])

xfit = xfit*np.max(rfu_data_mean)
yfit = yfit*np.max(cell_count)

# ax.plot(xfit,yfit)

ax.set_yscale('log')


# ax.annotate('$R^{2}$ = ' + str(round(res.rvalue**2,3)),(1,50000),fontsize=12)

# fig.savefig('cell_count_lin_regress.png',bbox_inches='tight')
# %% RFU to cell count

# res = stats.linregress(ydata,dilution_log[1:])

def rfu_to_cell_count(rfu):
    dilution_factor = 8.66 + (-1.812*10**-5)*rfu
    dilution_factor=2**(-dilution_factor)
    return dilution_factor*90000
# %%

# Graphs for fluorescence values vs. dilution using area under the curve from 0-60 min
# print('Graphs for fluorescence values vs. dilution using area under the curve from 0-60 min\n')
# fig3,ax = plt.subplots()

# col_indx=0

# curve_area_no_drug_avg = np.zeros(9)
# curve_area_no_drug = np.zeros((3,9))
# curve_area_drug = np.zeros((3, 9))

# for col in col_list[1:]:
#     ts_avg = np.zeros(len(time_vect))
    
#     row_indx = 0
#     for row in ['B','C','D']:
#         key = row + str(col)
#         ts = np.array(p_ab.data[key]).astype('float64') - bg
#         ts_spline = interp.InterpolatedUnivariateSpline(time_vect[:-1:]/60, ts[:-1:])
#         curve_area_no_drug[row_indx][col_indx] = ts_spline.integral(0, 60)
#         ts_avg += ts
#         row_indx += 1
#     ts_avg = ts_avg/3
#     ts_spline = interp.InterpolatedUnivariateSpline(time_vect[:-1:]/60, ts_avg[:-1:])

#     curve_area_no_drug_avg[col_indx] = ts_spline.integral(0, 60)

#     row_indx = 0
#     for row in ['E','F','G']:

#         key = row + str(col)
#         ts = np.array(p_ab.data[key]).astype('float64') - bg
#         ts_drug_spline = interp.InterpolatedUnivariateSpline(time_vect[:-1:]/60, ts[:-1:])
#         #ax_list3[row_indx].plot(xfine,ts_drug_spline(xfine, 0),color=cmap(col_indx/8),label=dilution[col_indx])
#         curve_area_drug[row_indx][col_indx] = ts_drug_spline.integral(0, 60)
#         row_indx+=1

#     col_indx+=1

# ax.plot(dilution_log[1:], curve_area_no_drug_avg, marker='o', markersize=3,color=cmap(0),label='no drug')
# ax.plot(dilution_log[1:], curve_area_drug[0], marker='o', markersize=3,color=cmap(0.33),label='1xMIC')
# ax.plot(dilution_log[1:], curve_area_drug[1], marker='o', markersize=3,color=cmap(0.66),label='10xMIC')
# ax.plot(dilution_log[1:], curve_area_drug[2], marker='o', markersize=3,color=cmap(0.99),label='100xMIC')

# ax.legend(frameon=False,fontsize=12)

# dilution_xlabel = ['$2^{' + str(d) + '}$' for d in dilution_log]
# ax.set_xticks(np.arange(1,10))
# ax.set_xticklabels(dilution_xlabel[1:])

# ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
# ax.tick_params(axis='both',labelsize=12)

# ax.set_ylabel('AUC RFU (a.u.)',fontsize=14)
# ax.set_xlabel('Dilution',fontsize=14)

# fig3.tight_layout()
# # %%
# all_auc_data = np.concatenate((curve_area_no_drug,curve_area_drug))
# ydata = np.mean(np.concatenate((curve_area_no_drug,curve_area_drug)),axis=0)
# yerr = np.std(np.concatenate((curve_area_no_drug,curve_area_drug)),axis=0)/(np.sqrt(6))

# res = stats.linregress(ydata,dilution_log[1:])

# def rfu_to_cell_count(rfu):
#     dilution_factor = 8.66 + (-1.812*10**-5)*rfu
#     dilution_factor=2**(-dilution_factor)
#     return dilution_factor*90000
# # %%
# K = (432/50)*10**5 # carrying capacity cells per uL
# cell_count = [K/(2**d) for d in dilution_log]
# cell_count_log = np.log(cell_count)

# # def logistic_eqn(x,y0,x_50,n):
# #     return y0/(1+(x/x_50)**n)

# def inv_logistic(y,y0,k,x50):
#     return x50/((y0/y - 1)**k)

# def inv_logistic_log(y,y0,k,x50):
#     return np.log10(x50/((y0/y - 1)**k))


# fig,ax = plt.subplots()

# norm_factor = np.max(ydata)
# auc_norm = ydata/np.max(ydata)
# yerr_norm = yerr/norm_factor
# cell_count_norm = cell_count/np.max(cell_count)

# ax.errorbar(auc_norm,cell_count_norm[1:],xerr=yerr_norm,yerr=None,fmt='o',color='black')
# # ax.set_yscale('log')

# # res = np.polyfit(auc_norm,cell_count_norm[1:],deg=3)

# xfit = np.arange(np.min(auc_norm),np.max(auc_norm),0.01)

# p0 = [1.1,1,0.01]
# bounds = [[1.001,0,0],[10,10,0.06]]
# popt,pcov = sciopt.curve_fit(inv_logistic,auc_norm,cell_count_norm[1:],maxfev=10000,p0=p0,bounds=bounds)
# # yfit = inv_logistic(xfit,popt[0],popt[1],popt[2])
# # yfit = inv_logistic(xfit,popt[0],1,0.05)
# yfit = inv_logistic(xfit,popt[0],popt[1],popt[2])
# # sfit = interp.InterpolatedUnivariateSpline(auc_norm,cell_count_norm[1:])



# # yfit = res[0]*xfit**2 + res[1]*xfit + res[2]
# ax.plot(xfit,yfit)
# ax.grid(visible=True,axis='both',which='both')
# # for data in all_auc_data:
# #     ax.plot(dilution_log[1:],data)

# # norm_factor = np.max(ydata)
# # ydata_norm = ydata/np.max(ydata)


# # popt,pcov = sciopt.curve_fit(logistic_eqn,cell_count_log[1:],ydata_norm)

# # cell_count_fit = np.arange(cell_count_log[1],cell_count_log[-1],0.1)
# # yfit = norm_factor*logistic_eqn(cell_count_fit,popt[0],popt[1],popt[2])

# # ax.plot(cell_count_fit,yfit,color='red')
# # # dilution_xlabel = ['$2^{' + str(d) + '}$' for d in dilution_log]

# # ax.set_xticks(np.arange(1,10))
# # ax.set_xticklabels(cell_count[1:])

# # ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
# # ax.tick_params(axis='both',labelsize=12)

# # ax.set_ylabel('AUC$_{60}$ (a.u.)',fontsize=14)
# # ax.set_xlabel('Dilution',fontsize=14)
# # %%
