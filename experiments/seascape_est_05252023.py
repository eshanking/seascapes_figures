#%%
import re
from fears.utils.AutoRate import Plate
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl

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

rate_est_lib = {}

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
fig,ax_list = plt.subplots(nrows=4,ncols=4,figsize=(15,12))
ax_list = ax_list.flatten()

for i in range(len(timeseries_dicts)):
    ts_dict = timeseries_dicts[i]
    ax = ax_list[i]
    for i in range(len(col_list)):
        ts_t = np.zeros((len(ts_dict['Time']),8))
        col = col_list[i]
        for j in range(len(row_list)):
            row = row_list[j]
            key = row+col
            if key in ts_dict.keys():
                ts_t[:,j] = ts_dict[key]
        ts_avg = np.mean(ts_t,axis=1)
        ts_err = np.std(ts_t,axis=1)/np.sqrt(8)
        ax.errorbar(ts_dict['Time'],ts_avg,yerr=ts_err,label=col,color=cmap(i/11))
    # ax.legend(frameon=False)

    # for row in row_list:
    #     col_indx = 0
    #     for col in col_list:
    #         key = row+col
    #         if key in ts_dict.keys():
    #             ax.plot(ts_dict['Time'],ts_dict[key],label=key,color=cmap(col_indx/12))
    #             col_indx += 1

#%%
    # get summary data

    # replicates = ['A','B','C','D','E','F','G','H']
    # conditions = np.linspace(1,12,num=12)
    # conditions = [str(int(c)) for c in conditions]

    # data_avg = {}
    # data_std = {}
    
    # gr_avg = []
    # gr_std = []

    # rate_est_dict = {}

    # for c in conditions:
    #     rate_est = []
    #     for r in replicates:

    #         if not (count == 6 and (r == 'A' or r == 'H')): 

    #             key = r+c

    #             # control_key = r+'12'

    #             # control_vect = np.array(timeseries_dict[control_key])

    #             od_vect = np.array(timeseries_dict[key])

    #             # od_vect = od_vect-control_vect

    #             od_vect= np.array(od_vect)

    #             timeseries_dict[key] = od_vect

    #             d,pcov = est_logistic_params(od_vect,t_vect,debug=False,mode='logistic',
    #                     normalize=False)

    #             logistic_params_dict[key] = d


    #             rate_est.append(d['gr'])
    #         else:
    #             key = r+c
    #             del timeseries_dict[key] # remove these two data points as outliers
    #         # r = rolling_regression(t_vect,od_vect)
    #         # rate_est.append(r)

    #     rate_est_dict[c] = rate_est

    #     # num_zeros = np.arghwere(rate_est)

    #     gr_avg.append(np.mean(rate_est))
    #     gr_std.append(np.std(rate_est))


    # grl_t = {'avg':gr_avg,
    #         'err':gr_std}
    
    # gr_lib[str(count)] = grl_t
    # rate_est_lib[str(count)] = rate_est_dict

    # all_timeseries[str(count)] = timeseries_dict
    # all_log_params[str(count)] = logistic_params_dict

    # count+=1