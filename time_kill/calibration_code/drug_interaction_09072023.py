from fears.utils import AutoRate
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp

plate_path = '../calibration_data/drug_interaction/EK_single_AB_constant_gain_20230907_114057.xlsx'
p_ab = AutoRate.Plate(plate_path,mode='single_measurement')
fluor_data = p_ab.od_data_to_dict(p_ab.data)

cols = np.arange(2,11)
rows = ['B','C','D','E','F','G']

drug_conc = [-3,-2,-1,0,1,2]
drug_conc = np.array(drug_conc)

fluor_mat = np.zeros((len(rows),len(cols)))

for row_indx,row in enumerate(rows):
    for col_indx,col in enumerate(cols):
        key = row + str(col)
        fluor_mat[row_indx,col_indx] = fluor_data[key]

fluor_avg = np.mean(fluor_mat,axis=1)
fluor_mat = fluor_mat/fluor_avg[0]
fluor_avg = fluor_avg/fluor_avg[0]
fluor_err = np.std(fluor_mat,axis=1)/np.sqrt(len(cols))

#%%
drug_conc_lin = [0,0.01,0.1,1,10,100]

fig,ax = plt.subplots()

ax.errorbar(drug_conc_lin,fluor_avg,yerr=fluor_err,fmt='o',capsize=5)
ax.set_xscale('symlog',linthresh=0.01)
# ax.set_xticks(drug_conc)
# drug_conc_plot = drug_conc = ['nd','$10^{-2}$','$10^{-1}$','$10^{0}$','$10^{1}$','$10^{2}$']
# ax.set_xticklabels(drug_conc_plot);

ax.set_xlabel('Drug Concentration (log$_{10}$)',fontsize=14)
ax.set_ylabel('Fluorescence (a.u.)',fontsize=14)
ax.tick_params(axis='both', which='major', labelsize=12)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
# %% spline fit

drug_conc_lin = np.array(drug_conc_lin)
drug_conc_lin = drug_conc_lin[1:]
fluor_avg = fluor_avg[1:]
fluor_err = fluor_err[1:]

drug_conc_log = np.log10(drug_conc_lin)

kind_list = ['linear', 'slinear', 'quadratic', 'cubic']

for kind in kind_list:
    f = interp.interp1d(drug_conc_log,fluor_avg,kind=kind)

    # drug_conc_lin = np.linspace(0.01,100,1000)
    drug_conc_fit = np.linspace(-2,2,1000)
    fluor_fit = f(drug_conc_fit)

    drug_conc_lin = 10**drug_conc_fit

    ax.plot(drug_conc_lin,fluor_fit,label=kind)

ax.set_ylim([1.01,1.04])
ax.set_xlim([0.01,100])
ax.legend(frameon=False,fontsize=10)

# %%

folder_path = '../experiment_data/tk_08302023'

