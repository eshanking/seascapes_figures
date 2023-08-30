#%%
from fears.utils import AutoRate
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy

plate_path = 'cell_count_calibration.xlsx'
p = AutoRate.Plate(plate_path)

carrying_cap = 10**6 # cells/uL

time_vect = np.array(p.data['Time [s]'])

row_list = ['B','C','D','E','F','G','H']
col_list = [2,3,4,5,6,7,8,9,10,11]
col_list = [str(c) for c in col_list]

row_pairs = [('B','C'),('D','E'),('F','G')]
data_avg = {}

cmap = mpl.colormaps['viridis']

condition_indx = 0

fig,ax = plt.subplots()

# for rp in row_pairs:
    # ax = ax_list[condition_indx]
dilution = []
for i in range(10):
    dilution.append(5**(-i))
rp = row_pairs[0]
col_indx = 0
for col in col_list:
    key1 = rp[0] + col
    key2 = rp[1] + col

    ts_avg = np.array(p.data[key1] + p.data[key2])/2
    ax.plot(time_vect/60,ts_avg,color=cmap(col_indx/10),label=str(round(np.log10(dilution[col_indx]),2)))

    data_avg[str(condition_indx)+col] = ts_avg

    col_indx += 1
ax.set_xlim(0,700)
ax.legend(frameon=False)
ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
ax.set_xlabel('Time (min)',fontsize=14)
ax.set_ylabel('RFU',fontsize=14)


#%% Time to 200000

thresh = 200000

condition = '0'

res = []
for col in col_list:
    key = condition + col
    ts = data_avg[key]
    indx = np.argwhere(ts>=thresh)[0][0]
    res.append(time_vect[indx]/60)

dilution = []
for i in range(10):
    dilution.append(5**(-i))

cell_count = np.array(dilution)*carrying_cap

fig,ax = plt.subplots()

ax.scatter(res,cell_count,marker='*')

ax.set_yscale('log')

ax.set_xlabel('Time to RFU 200000 (min)')
ax.set_ylabel('Cell count (uL$^{-1}$)');

cc_log = np.log10(cell_count[1:])
res_fit = np.array(res[1:])*60

lin_fit = np.polyfit(res_fit,cc_log,1)

time_to_thresh_fit = np.linspace(0,600*60)
cc_fit = lin_fit[0]*time_to_thresh_fit + lin_fit[1]

ax.plot(time_to_thresh_fit/60,10**cc_fit)

fig.savefig('time_to_2000_vs_cell_count.png',bbox_inches='tight')
# %% Fluorescence at 3600

condition = '0'
time_indx = np.argwhere(time_vect>=3600)[0][0]

res = []
for col in col_list:
    key = condition + col
    ts = data_avg[key]

    res.append(ts[time_indx])

fig,ax = plt.subplots()

ax.scatter(res[1:],cell_count[1:],marker='*')
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlabel('RFU')
ax.set_ylabel('Cell count (uL$^{-1}$)')


def fluor_v_dil_model(fluor,A,B,l):
    return A*np.log10(B*(fluor-l))

fluor_fit = np.log10(res[1:])
cc_fit = np.log10(cell_count[1:])

res_fit_censor = fluor_fit[cc_fit>2.5]
cc_fit_censor = cc_fit[cc_fit>2.5]

popt,pcov = scipy.optimize.curve_fit(fluor_v_dil_model,res_fit_censor,cc_fit_censor)

fluor_fit_plot = np.arange(np.min(res_fit_censor),np.max(res_fit_censor),0.01)
cc_fit_plot = fluor_v_dil_model(fluor_fit_plot,popt[0],popt[1],popt[2])

ax.plot(10**fluor_fit_plot,10**cc_fit_plot)

ax.annotate('log(N) = ' + str(round(popt[0],2)) + '*log[' + 
            str(round(popt[1],1)) + '*(log(RFU) - ' + str(round(popt[2],2)) + ')]', (10**4.7,10))


# dil_log = np.log10(dilution[1:])
# res_fit = np.log10(res[1:])

# lin_fit = np.polyfit(dil_log,res_fit,3)

# dil_fit = np.linspace(-7,-1)
# fluor_fit = lin_fit[0]*dil_fit**3 + lin_fit[1]*dil_fit**2 + lin_fit[2]*dil_fit + lin_fit[3]

# fluor_fit = np.linspace(20000,300000)
# dil_fit = lin_fit[0]*fluor_fit**2 + lin_fit[1]*fluor_fit + lin_fit[2]

# ax.plot(10**dil_fit,10**fluor_fit)

# %%
