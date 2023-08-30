#%%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.animation as animation
from fears.utils import AutoRate
import scipy
import os
import pandas as pd
import pickle

#%%
def get_polyfit_result(x,res):
    y = 0
    n = len(res)-1
    indx = 0
    for r in res:
        y += r*x**(n-indx)
        indx+=1
    return y

data_folder = 'calibration'

exp_files = os.listdir(data_folder)

exp_files = [f for f in exp_files if f[-5:] == '.xlsx']

exp_files.sort()

# p = AutoRate.Plate(data_file)

final_plate = np.zeros((8,10))

row_list = ['A','B','C','D','E','F','G','H']
col_list = [2,3,4,5,6,7,8,9,10,11]
col_list = [str(col) for col in col_list]

plate_list = []

for ef in exp_files:
    path_t = os.getcwd() + os.sep + data_folder + os.sep + ef
    p = AutoRate.Plate(path_t)   
    p.data = p.data.drop(p.data.index[0:2]) # remove first datapoint because of a systematic drop in fluorescence
    p.data = p.data.drop(p.data.index[-1]) # remove the last datapoint because of incomplete data

    dt = p.get_start_time()
    start_time = 60*((60*dt.hour) + dt.minute) + dt.second
    p.data['Time [s]'] = p.data['Time [s]'] + start_time

    plate_list.append(p)

#%% combine plates
data = plate_list[0].data

for p in plate_list[1:]:
    data = pd.concat((data,p.data))

data = data.sort_values(by='Time [s]')

data['Time [s]'] = data['Time [s]'] - np.min(data['Time [s]'])

#%%

for row_indx in range(8):
    for col in range(1,11):
        key = row_list[row_indx] + str(col+1)
        ts = np.array(data[key])

        final_plate[row_indx,col-1] = ts[-1]

# final_plate = final_plate - np.min(final_plate)
final_plate = final_plate/np.min(final_plate)

fig,ax = plt.subplots()
im = ax.imshow(final_plate)
fig.colorbar(im,ax=ax)
#%%

fig,ax = plt.subplots()

# dc = [1000,100,10,1,0.1,0.01,0.001,0]
dc = [0,0.001,0.01,0.1,1,10,100,1000]

bp = ax.boxplot(final_plate.T,labels=dc);
ax.set_ylabel('Normalized fluorescence',fontsize=14)
ax.set_xlabel('Drug concentration (ug/mL)',fontsize=14)

dc_log = [-4,-3,-2,-1,0,1,2,3]

mean_fluor = np.mean(final_plate,axis=1)

res = np.polyfit(dc_log,mean_fluor,3)

dc_fit = np.linspace(-4,3,num=100)

yfit = get_polyfit_result(dc_fit,res)

xt = ax.get_xticks()
dc_fit_plot = np.linspace(xt[0],xt[-1],num=100)

ax.plot(dc_fit_plot,yfit,color='red')
# %%
cmap = mpl.colormaps['viridis']
fig,ax_list = plt.subplots(nrows=10,figsize=(3,8))

time_vect = data['Time [s]']

for row_indx in range(8):
    for col in range(1,11):
        key = row_list[row_indx] + str(col+1)
        ts = np.array(data[key])
        ax_list[col-1].plot(time_vect,ts,color=cmap(row_indx/8))
        ax_list[col-1].set_xlim(0,3600)
# %%

spline_fit = scipy.interpolate.CubicSpline(dc_log,mean_fluor)

fig,ax = plt.subplots()
bp = ax.boxplot(final_plate.T,labels=dc);

xt = ax.get_xticks()
dc_fit_plot = np.linspace(xt[0],xt[-1],num=100)

ax.plot(dc_fit_plot,spline_fit(dc_fit),color='red',linewidth=2)

# %%

# def get_ctx_norm_factor(data,plate_dc,conc,t,debug=False):

#     fluor_data = np.zeros((8,10))

#     row_list = ['A','B','C','D','E','F','G','H']

#     # get index

#     time = np.array(data['Time [s]'])

#     indx = np.argwhere(time>t)[0][0]

#     for row_indx in range(8):
#         for col in range(1,11):
#             key = row_list[row_indx] + str(col+1)
#             ts = np.array(data[key])
#             fluor_data[row_indx,col-1] = ts[indx]

#     fluor_data = fluor_data/np.min(fluor_data)

#     mean_fluor = np.mean(fluor_data,axis=1)

#     spline_fit = scipy.interpolate.CubicSpline(plate_dc,mean_fluor)

#     if debug:
#         fig,ax = plt.subplots()
#         bp = ax.boxplot(fluor_data.T,labels=dc);

#         xt = ax.get_xticks()
#         dc_fit_plot = np.linspace(xt[0],xt[-1],num=100)

#         ax.plot(dc_fit_plot,spline_fit(dc_fit),color='red',linewidth=2)

#     return spline_fit(conc)

# %% animation

# fig, ax = plt.subplots()
# xdata, ydata = [], []
# ln, = ax.plot([], [], 'ro')

# def init():
#     ax.set_xlim(dc_fit_plot[0],dc_fit_plot[-1])
#     ax.set_ylim(1,2)
#     ax.set_xlabel('Log drug concentration')
#     ax.set_ylabel('Normalized fluorescence')
#     return ln,

# def update(frame):
#     fluor_data = np.zeros((8,10))
#     for row_indx in range(8):
#         for col in range(1,11):
#             key = row_list[row_indx] + str(col+1)
#             ts = np.array(data[key])
#             fluor_data[row_indx,col-1] = ts[frame]

#     fluor_data = fluor_data/np.min(fluor_data)

#     mean_fluor = np.mean(fluor_data,axis=1)

#     spline_fit = scipy.interpolate.CubicSpline(dc_log,mean_fluor)
#     ln.set_data(dc_fit_plot,spline_fit(dc_fit))
#     ax.set_title(str(frame))
#     ax.set_xticks(ax.get_xticks(),['ND',-3,-2,-1,0,1,2,3])
#     return ln,

# ani = animation.FuncAnimation(fig, update, frames=788,
#                     init_func=init, blit=True)

# writervideo = animation.FFMpegWriter(fps=60) 
# ani.save('spline_through_time.mov', writer=writervideo)
# %%

class Normalizer:
    def __init__(self,data_folder,plate_dc):
        
        self.plate_dc = plate_dc

        exp_files = os.listdir(data_folder)

        exp_files = [f for f in exp_files if f[-5:] == '.xlsx']

        exp_files.sort()

        self.exp_files = exp_files

        # p = AutoRate.Plate(data_file)

        self.row_list = ['A','B','C','D','E','F','G','H']
        col_list = [2,3,4,5,6,7,8,9,10,11]
        self.col_list = [str(col) for col in col_list]

        plate_list = []

        for ef in exp_files:
            path_t = os.getcwd() + os.sep + data_folder + os.sep + ef
            p = AutoRate.Plate(path_t)   
            p.data = p.data.drop(p.data.index[0:2]) # remove first datapoint because of a systematic drop in fluorescence
            p.data = p.data.drop(p.data.index[-1]) # remove the last datapoint because of incomplete data
            dt = p.get_start_time()
            start_time = 60*((60*dt.hour) + dt.minute) + dt.second
            p.data['Time [s]'] = p.data['Time [s]'] + start_time

            plate_list.append(p)

        # combine plates
        data = plate_list[0].data

        for p in plate_list[1:]:
            data = pd.concat((data,p.data))

        data = data.sort_values(by='Time [s]')

        data['Time [s]'] = data['Time [s]'] - np.min(data['Time [s]'])

        self.data = data

    def get_ctx_norm_factor(self,conc,t,debug=False):
        """Returns a scalar to normalize fluorescence values to compare to the 
           no drug condition.

        Args:
            conc (float): log of concentration (concentration must be > 0)
            t (time): time of measurements (seconds)
            debug (bool, optional): If True, creates a box plot and 
            corresponding spline fit of the calibration data. Defaults to False.

        Returns:
            float: normalization factor
        """

        fluor_data = np.zeros((8,10))

        row_list = ['A','B','C','D','E','F','G','H']

        # get index

        time_vect = np.array(self.data['Time [s]'])

        indx = np.argwhere(time_vect>=t)[0][0]

        for row_indx in range(8):
            for col in range(1,11):
                key = row_list[row_indx] + str(col+1)
                ts = np.array(self.data[key])
                fluor_data[row_indx,col-1] = ts[indx]

        fluor_data = fluor_data/np.min(fluor_data)

        mean_fluor = np.mean(fluor_data,axis=1)

        spline_fit = scipy.interpolate.CubicSpline(self.plate_dc,mean_fluor)

        if debug:
            fig,ax = plt.subplots()
            bp = ax.boxplot(fluor_data.T,labels=dc);

            xt = ax.get_xticks()
            dc_fit_plot = np.linspace(xt[0],xt[-1],num=100)

            ax.plot(dc_fit_plot,spline_fit(dc_fit),color='red',linewidth=2)

        return float(spline_fit(conc))

# %%

norm = Normalizer('calibration',dc_log)
filehandler = open('normalizer.p', 'wb') 
pickle.dump(norm, filehandler)

# %%
