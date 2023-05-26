#%%
import re
from fears.utils.AutoRate import Plate
from fears.utils import plotter
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.optimize as sciopt
import scipy.stats as stats

# reminder: col 12 is the background condition

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

def get_plate_paths(folder_path):
    """Gets plate data paths
    Returns:
        list: list of plate data paths
    """
    plate_files = os.listdir(path=folder_path)

    #Need to make sure we are only attempting to load .csv or .xlsx data
    plate_files = [i for i in plate_files]

    plate_files.sort()

    plate_data_paths = []

    for pf in plate_files:
        if pf != '.DS_Store':
            plate_path = folder_path + os.sep + pf
            plate_data_paths.append(plate_path)

    plate_data_paths.sort(key=natural_keys)
    return plate_data_paths

def get_data_file_paths(plate_path):
    files = os.listdir(path=plate_path)

    #Need to make sure we are only attempting to load .csv or .xlsx data
    files = [i for i in files if ('.csv' in i) or ('.xlsx' in i)]

    files.sort()

    file_data_paths = []

    for pf in files:
        if pf != '.DS_Store':
            file_path = plate_path + os.sep + pf
            file_data_paths.append(file_path)

    return file_data_paths

#%%
row_list = ['A','B','C','D','E','F','G','H']
col_list = ['1','2','3','4','5','6','7','8','9','10','11','12']
drug_conc = [10000,2000,400,80,16,3.2,0.64,0.128,0.0256,0.00512,0,'control']
folder_path = 'data/08312022'

plate_paths = get_plate_paths(folder_path)

def get_data_file_paths(plate_path):
    files = os.listdir(path=plate_path)

    #Need to make sure we are only attempting to load .csv or .xlsx data
    files = [i for i in files if ('.csv' in i) or ('.xlsx' in i)]

    files.sort()

    file_data_paths = []

    for pf in files:
        if pf != '.DS_Store':
            file_path = plate_path + os.sep + pf
            file_data_paths.append(file_path)

    return file_data_paths

#%%
# fig,ax_list = plt.subplots(ncols=4,nrows=4,figsize=(15,12))

bg_col = '12'

count = 0

gr_lib = {}

# all_timeseries = {}
# all_log_params = {}
# list of dicts representing each plate. Keys of the dicts are individual wells
timeseries_dicts = []

for pp in plate_paths:
# for pp in [plate_paths[9]]:

    row = int(np.floor(count/4))
    col = int(np.mod(count,4))

    # ax = ax_list[row,col]

    # fig,ax = plt.subplots()

    data_paths0 = get_data_file_paths(pp)

    timeseries_dict = {}
    timeseries_dict['Time'] = []
    logistic_params_dict = {}

    for p in data_paths0:

        plate = Plate(p,mode='single_measurement')
        data = plate.data
        data_dict = plate.od_data_to_dict(data)

        # df = pd.read_excel(p)
        t = plate.get_start_time()

        data_dict['Time'] = t

        # estimate background
        bg_est = 0
        for row in row_list:
            key = row + bg_col
            bg_est += data_dict[key]
        bg_est = bg_est/len(row_list)

        # bg = get_background(data_dict,bg_keys)
        for key in data_dict:
            if key != 'Time':
                if key in timeseries_dict.keys():
                    od = data_dict[key] - bg_est
                    timeseries_dict[key].append(od)
                else:
                    od = data_dict[key] - bg_est
                    timeseries_dict[key] = [od]
        timeseries_dict['Time'].append(t)

    # sort out time
    t_vect = timeseries_dict['Time']
    t0 = t_vect[0]
    t_vect = [(t-t0).total_seconds() for t in t_vect]
    timeseries_dict['Time'] = t_vect
    timeseries_dicts.append(timeseries_dict)

#%% visualize od data
cmap = mpl.cm.get_cmap('viridis')
# fig,ax_list = plt.subplots(nrows=4,ncols=4,figsize=(15,12))

auc_seascape = {}
max_auc = 0

for i in range(len(timeseries_dicts)):
# for i in range(1):
    # fig,ax_list = plt.subplots(nrows=3,ncols=4,sharex=True,sharey=True,figsize=(10,8))
    # ax_list = ax_list.flatten()
    fig,ax = plt.subplots()
    ts_dict = timeseries_dicts[i]
    time = np.array(ts_dict['Time'])/60
    
    auc_v_dc = []
    auc_v_dc_err = []
    for j in range(len(col_list)):
        # ax = ax_list[j]
        ts_t = np.zeros((len(ts_dict['Time']),8))
        col = col_list[j]
        auc_list = []
        for k in range(len(row_list)):
            row = row_list[k]
            key = row+col
            if key in ts_dict.keys():
                ts_t[:,k] = ts_dict[key]
                # auc_list.append(np.trapz(ts_t[:,k],ts_dict['Time']))
                auc_list.append(np.trapz(ts_t[:,k],time))
                # ax.plot(ts_dict['Time'],ts_t[:,k])
        # ts_t = ts_t - np.min(ts_t) + 10**-10
        ts_avg = np.mean(ts_t,axis=1)
        # ts_log = np.log10(ts_t)
        # ts_log_avg = np.mean(ts_log,axis=1)
        # ts_log_err = np.std(ts_log,axis=1)/np.sqrt(8)
        ts_err = np.std(ts_t,axis=1)/np.sqrt(8)

        # ax.errorbar(ts_dict['Time'],ts_avg,yerr=ts_err,label=col,color=cmap(i/11))
        # calculate AUC with ts_t
        auc_mean = np.mean(auc_list)
        auc_err = np.std(auc_list)/np.sqrt(8)
        # auc_seascape_t = {'avg':np.mean(auc_list),'err':np.std(auc_list)/np.sqrt(8)}
        # auc_seascape[col] = auc_seascape_t

        auc_v_dc.append(auc_mean)
        auc_v_dc_err.append(auc_err)

        if auc_mean > max_auc:
            max_auc = auc_mean
    
    ax.errorbar(drug_conc[0:-1],auc_v_dc[0:-1],yerr=auc_v_dc_err[0:-1],fmt='o')
    # fig.tight_layout()
    fig.suptitle('genotype = ' + str(i),y=1.05)
    ax.set_xscale('log')
    
    auc_seascape_t = {'avg':auc_v_dc,'err':auc_v_dc_err}
    auc_seascape[i] = auc_seascape_t

#%% fit pharmacodynamic model

def hill_fn(conc,gmax, gmin, hc, ic_50):
    
    y = gmax + ((gmin - gmax) * conc**hc) / (ic_50**hc + conc**hc)
    return y

seascape_lib = {}
gmax_list = []
g0_list = []
ic_50_list = []

fig,ax_list = plt.subplots(ncols=4,nrows=4,figsize=(12,10),sharex=True,sharey=True)
ax_list = ax_list.flatten()

normalized_auc_seascape = {}

for key in auc_seascape:
    ax = ax_list[key]    
    auc_seascape_t = auc_seascape[key]
    auc_v_dc = np.array(auc_seascape_t['avg'])
    auc_v_dc_err = np.array(auc_seascape_t['err'])

    # normalize by max auc
    auc_v_dc = auc_v_dc/max_auc
    auc_v_dc_err = auc_v_dc_err/max_auc

    normalized_auc_seascape_t = {'avg':auc_v_dc,'err':auc_v_dc_err}
    normalized_auc_seascape[key] = normalized_auc_seascape_t

    # estimate gmax and gmin
    gmax = np.max(auc_v_dc[0:-1])
    gmin = np.min(auc_v_dc[0:-1])

    p0 = [gmax, gmin, 1, 1]
    bounds = ([gmin,-1,0,0],[2,gmax,np.inf,np.inf])

    popt, pcov = sciopt.curve_fit(hill_fn, drug_conc[0:-1], auc_v_dc[0:-1], p0=p0)
    gmax = popt[0]
    gmin = popt[1]
    hc = popt[2]
    ic_50 = popt[3]
    rate_est = {'gmax':gmax,'gmin':gmin,'hc':hc,'ic_50':ic_50}
    seascape_lib[key] = rate_est
    xfit = np.logspace(-3,4,100)
    yfit = hill_fn(xfit,gmax,gmin,hc,ic_50)
    ax.errorbar(drug_conc[0:-1],auc_v_dc[0:-1],yerr=auc_v_dc_err[0:-1],fmt='o',markersize=3)
    ax.plot(xfit,yfit)
    # fig.suptitle('genotype = ' + str(key))
    ax.set_xscale('log')

    gmax_list.append(gmax)
    g0_list.append(auc_v_dc[-2])
    ic_50_list.append(ic_50)

# %%

fig,ax = plt.subplots(figsize=(5,5))

def int_to_binary(num,n_genotype=16):
    """
    Converts an integer to binary representation with the number of 
    digits equal to the number of alleles in the model.

    Parameters
    ----------
    num : int
        Number to be converted.

    Returns
    -------
    str
        Binary representation.

    """
    pad = int(np.log2(n_genotype))
    return bin(num)[2:].zfill(pad)

cmap = mpl.cm.get_cmap('magma',6)
cmap = cmap.colors
indx = 0

mut_list = []

for i in range(16):

    ic50 = np.log10(ic_50_list[i])
    g_drugless = g0_list[i]
    key_bin = int_to_binary(int(i))

    num = 0
    for s in key_bin:
        num+=int(s)

    mut_list.append(num)
    ax.scatter(ic50,g_drugless,marker='o',s=400,facecolor=cmap[num+1],
                edgecolors='w',label=int(num))
    # ax4.annotate(key,(ic50-0.15,g_drugless-0.001),fontsize=12)
    ax.annotate(i,(ic50,g_drugless),fontsize=12,ha='center',va='center')

# ax.set_ylim(0.06,0.115)
ax.set_xlim(-3,4)
ax.set_ylabel('Normalized AUC',fontsize=14)
ax.set_xlabel('Log IC50 (ug/mL)',fontsize=14)
ax.tick_params(axis='both', labelsize=13)

handles, labels = ax.get_legend_handles_labels()

unique_labels = sorted(set(labels))
labels = np.array(labels)
unique_handles = []

for lab in unique_labels:
    indx = np.argwhere(labels==lab)
    indx = indx[0][0]
    unique_handles.append(handles[indx])

ax.legend(unique_handles,unique_labels,loc = (1,0),frameon=False,
             fontsize=12,title='$n_{mut}$')


# %%
fig8, ax_list = plt.subplots(ncols=2,figsize=(8,4))

ax = ax_list[0]

# wt_gr = seascape_lib['0']['g_drugless']
gr_list_norm= g0_list/g0_list[0]

ax.scatter(mut_list,gr_list_norm)

ax.set_xticks([0,1,2,3,4])

res = stats.linregress(mut_list,gr_list_norm)

x = np.arange(5)
y = res.slope*x + res.intercept

ax.plot(x,y,color='orange',linewidth=2)

ax.set_xlabel('# mutations',fontsize=14)
ax.set_ylabel('Normalized growth rate',fontsize=14)
ax.tick_params(axis='both', labelsize=13)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

ann = '$r^2$ = ' + str(round(res.rvalue**2,2)) + '\n$p$ = ' + str(round(res.pvalue,3))
ax.annotate(ann,(2.5,1.2),fontsize=12)

# Mutations vs IC50
ax = ax_list[1]

# wt_ic50 = seascape_lib['0']['ic50']

# ic50_list_norm = [10**x for x in ic50_list]

ic50_list_norm = ic_50_list/(ic_50_list[0])

ic50_list_norm = [np.log10(x) for x in ic50_list_norm]

ax.scatter(mut_list,ic50_list_norm)

res = stats.linregress(mut_list,ic50_list_norm)

x = np.arange(5)
y = res.slope*x + res.intercept

ax.plot(x,y,color='orange',linewidth=2)

ax.set_xticks([0,1,2,3,4])

ax.set_xlabel('# mutations',fontsize=14)
ax.set_ylabel('Normalized IC50',fontsize=14)
ax.tick_params(axis='both', labelsize=13)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

ann = '$r^2$ = ' + str(round(res.rvalue**2,2)) + '\n$p$ = ' + str(round(res.pvalue,3))
ax.annotate(ann,(2.5,0.2),fontsize=12)

fig8.tight_layout()
# %% plot seascape

fig,ax = plt.subplots(figsize=(10,6))

cc = plotter.gen_color_cycler().by_key()
# ax.set_prop_cycle(cc)

drug_conc_plot = ['10^$' + str(round(np.log10(x),1)) + '$' for x in drug_conc[0:-2]]

drug_conc_plot.append('nd')
x = np.linspace(-4,4)

for key in seascape_lib:
    color = cc['color'][key]
    linestyle = cc['linestyle'][key]
    rate_est = seascape_lib[key]
    gmax = rate_est['gmax']
    gmin = rate_est['gmin']
    hc = rate_est['hc']
    ic_50 = rate_est['ic_50']
    xfit = np.logspace(-3,4,100)
    xfit = np.concatenate(([0],xfit))
    yfit = hill_fn(xfit,gmax,gmin,hc,ic_50)

    ax.plot(xfit,yfit,linestyle,label=key,color=color,linewidth=2)

    auc = normalized_auc_seascape[key]['avg']
    auc_err = normalized_auc_seascape[key]['err']
    # ax.errorbar(drug_conc[0:-1],auc[0:-1],yerr=auc_err[0:-1],
    #             fmt='o',markersize=3,color=color)

ax.set_xscale('symlog',linthresh=10**-3)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xlabel('Drug concentration (ug/mL)',fontsize=14)
ax.set_ylabel('Normalized AUC',fontsize=14)
ax.tick_params(axis='both', labelsize=12)
ax.legend(title='genotype',fontsize=10,frameon=False,ncol=1,loc=(1,0.1))

# %%
mic_list = [0.088,1.4,0.063,32,0.13,3.6*10**2,0.18,3.6*10**2,0.088,23,1.4,
            3.6*10**2,1.4,2.1*10**3,0.8,2.9*10**3]