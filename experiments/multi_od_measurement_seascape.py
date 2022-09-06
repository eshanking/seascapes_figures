import os
import pandas as pd
import datetime
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as colors
import matplotlib.cm as mplcm
import scipy.optimize as sciopt

def logistic_growth_curve(t,r,p0,k):
    """Logistic growth equation

    Args:
        t (float): time
        r (float): growth rate
        p0 (float): starting population size
        k (float): carrying capacity

    Returns:
        float: population size at time t
    """
    p = k/(1+((k-p0)/p0)*np.exp(-r*t))

    return p

def est_logistic_params(growth_curve,t,debug=False,sigma=None):
    """Estimates growth rate from OD growth curve

    Args:
        growth_curve (list or numpy array): vector of OD data
        t (list or numpy array, optional): Time vector. If None, algorithm assumes each time step is 1 s. Defaults to None.

    Returns:
        dict: Dict of logistic growth curve paramters
    """

    p0 = [10**-3,0.1,0.5] # starting parameters

    popt, pcov = sciopt.curve_fit(logistic_growth_curve,
                                        t,growth_curve,p0=p0,sigma=sigma)
    
    rate_indx = 0 # index of the optimized data points referring to the growth rate
    p0_indx = 1 # index of the optimized data points referring to initial population size
    carrying_cap_indx = 2 # index of the optmized data points referring to carrying capacity

    r = popt[rate_indx]
    p0 = popt[p0_indx]
    cc = popt[carrying_cap_indx]

    min_carrying_cap = 0.2

    if r < 0: # if the growth rate is negative
        r = 0
    # if cc < p0: # if the carrying capacity is less than the initial population size
    #     r = 0
    # if cc < min_carrying_cap: # if the carrying cap is less the minimum threshold
    #     r = 0

    d = {'gr':r,
            'OD_0':p0,
            'OD_max':cc}   

    if debug:
        fig,ax = plt.subplots()

        ax.plot(t,growth_curve)

        est = logistic_growth_curve(t,popt[0],popt[1],popt[2])
        
        ax.plot(t,est)
        # print(popt[0])
        p0 = round(popt[1]*10**5)/10**5
        k = round(popt[2]*10**5)/10**5
        r = round(r*10**5)/10**5
        title = 'rate = ' + str(r*(60**2)) + ' cc = ' + str(k)
        ax.set_title(title)   

    return d,pcov

def get_background(data_dict,bg_keys):
    avg = 0
    for key in bg_keys:
        avg+=data_dict[key]
    avg = avg/len(bg_keys)
    return avg

def od_data_to_dict(df):
    """Takes an OD data file and returns a dict with each data well as a key

    OD data files refers to the type of data that is a single OD reading in time (i.e
    not timeseries data).

    Args:
        df (pandas dataframe): Parsed OD data file

    Returns:
        dict: Dict of key-value pairs where each key is a well location and each value
        is the OD reading.
    """
    rownum = 0
    d = {}

    for k0 in df['Rows']:
        for k1 in df.columns[1:]:
            if k0.isnumeric():
                key = k1+k0
            else:
                key = k0+k1
            d[key] = df[k1][rownum]
        rownum+=1

    return d                    

def parse_od_data_file(df):
    """Loads the raw OD data files and strips metadata

    OD data files refers to the type of data that is a single OD reading in time (i.e
    not timeseries data).

    Args:
        data_path (str): path to data file

    Returns:
        pandas dataframe: Formatted dataframe with metadata stripped
    """

    # get the first column as an array
    col_0 = df.columns[0]
    col_0_array = np.array(df[col_0])

    data_start_indx = np.argwhere(col_0_array=='<>')
    data_start_indx = data_start_indx[0][0]

    # the data end index should be the first NaN after the data start index
    col_0_array_bool = [pd.isna(x) for x in col_0_array]
    data_end_indx = np.argwhere(col_0_array_bool[data_start_indx:])
    data_end_indx = data_end_indx[0][0] + data_start_indx - 1

    df_filt = df.loc[data_start_indx:data_end_indx,:]

    # fix the columns
    i = 0
    columns = list(df_filt.iloc[0])
    columns_t = []
    for c in columns:
        if type(c) is not str:
            columns_t.append(str(int(c)))
        else:
            columns_t.append(c)
        i+=1
    
    columns_t[0] = 'Rows'

    df_filt.columns = columns_t

    df_filt = df_filt.drop(df_filt.index[0])
    df_filt = df_filt.reset_index(drop=True)

    return df_filt

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

def get_start_time(df,col=4):

    # first start time is shaking, so we take the second (start of scan)
    f = df[df == 'Start Time'].stack().index.tolist()[1]

    row = f[0]
    date_time = df.iloc[row,col]

    yr = int(date_time[0:4])
    mon = int(date_time[5:7])
    day = int(date_time[8:10])

    hr = int(date_time[11:13])
    min = int(date_time[14:16])
    sec = int(date_time[17:19])

    dt = datetime.datetime(yr,mon,day,hour=hr,minute=min,second=sec)


    return dt

num_colors = 12
cm = cm.get_cmap('viridis',12)

cNorm  = colors.Normalize(vmin=0, vmax=num_colors-1)
scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)

bg_keys = ['A12','B12','C12','D12','E12','F12','G12','H12']
drug_conc = [10000,2000,400,80,16,3.2,0.64,0.128,0.0256,0.00512,0,0]
folder_path = '/Users/eshanking/repos/seascapes_figures/data/08312022'

plate_paths = get_plate_paths(folder_path)

fig,ax_list = plt.subplots(ncols=4,nrows=4,figsize=(10,8))

count = 0

for pp in plate_paths:

    row = int(np.floor(count/4))
    col = int(np.mod(count,4))

    ax = ax_list[row,col]

    # fig,ax = plt.subplots()

    data_paths0 = get_data_file_paths(pp)

    timeseries_dict = {}
    timeseries_dict['Time'] = []

    for p in data_paths0:

        df = pd.read_excel(p)
        t = get_start_time(df)
        df = parse_od_data_file(df)

        data_dict = od_data_to_dict(df)
        data_dict['Time'] = t

        # bg = get_background(data_dict,bg_keys)
        for key in data_dict:
            if key != 'Time':
                if key in timeseries_dict.keys():
                    od = data_dict[key]
                    timeseries_dict[key].append(data_dict[key])
                else:
                    od = data_dict[key]
                    timeseries_dict[key] = [od]
        timeseries_dict['Time'].append(t)

    # sort out time
    t_vect = timeseries_dict['Time']
    t0= t_vect[0]
    t_vect = [(t-t0).total_seconds() for t in t_vect]

    # get summary data

    replicates = ['A','B','C','D','E','F','G','H']
    conditions = np.linspace(1,12,num=12)
    conditions = [str(int(c)) for c in conditions]

    data_avg = {}
    data_std = {}

    for t in range(len(data_paths0)):
        for c in conditions:
            d = []
            for r in replicates:
                key = r + c
                od = timeseries_dict[key][t]
                d.append(od)

            if c in data_avg.keys():
                data_avg[c].append(np.mean(d))
                data_std[c].append(np.std(d))
            else:
                data_avg[c] = [np.mean(d)]
                data_std[c] = [np.std(d)]

    # for key in timeseries_dict:
    #     if key != 'Time':
    #         drug_indx = int(key[1:])
    #         dc = drug_conc[drug_indx-1]
    #         color = cmap[drug_indx-1]
    #         ax.scatter(t_vect,timeseries_dict[key],color=color,label=str(dc))

    for c in conditions:

        # curve fit

        d,pcov = est_logistic_params(data_avg[c],t_vect)

        drug_indx = int(c)-1
        dc = drug_conc[drug_indx]
        color = scalarMap.to_rgba([drug_indx])
        ax.errorbar(t_vect,data_avg[c],yerr=data_std[c],label=c,color=color,fmt='*')

        # plot curve fit
        # if d['gr'] != 0:
        cf = [logistic_growth_curve(t,d['gr'],d['OD_0'],d['OD_max']) for t in t_vect]
        ax.plot(t_vect,cf,color=color)
        # d_t = d


    handles, labels = ax.get_legend_handles_labels()

    # ax.legend(handles[0:12], labels[0:12])
    ax.set_title(str(count))
    count += 1

plt.tight_layout()
fig.savefig('logistic_growth_fit.pdf')