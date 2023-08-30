#%%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.integrate import odeint
import scipy
from fears.utils import AutoRate

data_file = 'EK_AB_20230308_125147.xlsx'

ab_ctx_cal_file = ''

row_list = ['A','B','C','D','E','F','G','H']

def get_cells_from_od(od):
    return (1.48*10**7)*od**2 + (1.62*10**7)*od -4.5*10**5

def rolling_average(x,N):
    
    indx = 0
    res = []
    while indx+N < len(x):
        xt = x[indx:indx+N]
        res.append(np.nanmean(xt))
        indx+=1
    
    for i in range(len(x[indx:])):
        res.append(np.nanmean(x[indx+i:]))

    res = np.array(res)
    x = np.array(x)
    res[x == 0] = 0

    return res

def kinetic_eqn(y,t,ka,kd,N):
    # print(y)
    r,R = y
    ka = 10**ka
    kd = 10**kd
    N = 10**N
    dydt = [-ka*r*N,ka*r*N - kd*R]
    return dydt

def kinetic_sol(t,y0,ka,kd,N):
    y = odeint(kinetic_eqn,y0,t,args=(ka,kd,N))
    return y

# Define a global error function that sums up the squared errors for each data set 
def global_loss(params,ydata,t):

    ka,kd,Ns = params[0],params[1],params[2:]
    indx = 0
    loss = 0

    for yd in ydata:
        y_t = kinetic_sol(t,[1,0],ka,kd,Ns[indx])[:,1]
        loss += (np.array(y_t - yd)**2).sum()
        indx+=1
    
    return loss

def global_loss_known_N(params,ydata,t,Ns):

    ka,kd = params[0],params[1]
    indx = 0
    loss = 0

    for yd in ydata:
        y_t = kinetic_sol(t,[1,0],ka,kd,Ns[indx])[:,1]
        loss += (np.array(y_t - yd)**2).sum()
        indx+=1
    
    return loss

def get_polyfit_result(x,res):
    y = 0
    n = len(res)-1
    indx = 0
    for r in res:
        y += r*x**(n-indx)
        indx+=1
    return y

p = AutoRate.Plate(data_file)

# fig,ax_list = plt.subplots(nrows=11,figsize=(3,8))

# # ax_list = ax_list.reshape(-1)

time = p.data['Time [s]']
cmap = mpl.colormaps['viridis']

# for col in range(1,12):
#     ax = ax_list[col-1]
#     row_indx = 0
#     for row in row_list:
#         key = row + str(col)
#         ax.plot(time,p.data[key],color=cmap(row_indx/8),linewidth=2)
#         row_indx+=1
#%% Drug vs ctx calibration

fig,ax = plt.subplots()

col = '12'
row_indx = 0
bg_data = p.data['A12']
dc = [0,0.1,1,10,100,1000,10000]

fluor_v_dc = []

for row in row_list[0:-1]:
    key = row+col
    ts = p.data[key]
    ts = np.array(np.divide(ts,bg_data))
    ax.plot(time,ts,color=cmap(row_indx/7),linewidth=2,label=dc[row_indx])
    fluor_v_dc.append(ts[-1])
    row_indx+=1

ax.legend()

#%%

fig,ax = plt.subplots()

dc_log = [-2,-1,0,1,2,3,4]

res_f_v_dc = np.polyfit(dc_log,fluor_v_dc,4)

dc_fit = np.linspace(-2,4,num=100)
# fluor_fit = res[0]*dc_fit**4 + res[1]*dc_fit**3 + res[2]*dc_fit**2 + res[3]*dc_fit + res[4]
fluor_fit = get_polyfit_result(dc_fit,res_f_v_dc)

ax.plot(dc_log,fluor_v_dc,'*',markersize=10)
ax.plot(dc_fit,fluor_fit,color='red',linewidth=2)
ax.set_ylabel('Relative fluorescence',fontsize=14)
ax.set_xlabel('log drug concentration',fontsize=14)
# ax.set_xscale('log')
# ax.set_xlim(10000,30000)
#%%
fig,ax = plt.subplots()

row_indx = 0
col = '11'

max_val = np.max(p.data['A'+col])

for row in row_list[0:-1]:
    key = row+col
    dc = dc_log[row_indx]
    norm_coeff = get_polyfit_result(dc,res_f_v_dc)
    ts = np.array(p.data[key])
    ts = ts/norm_coeff
    ts = ts/max_val
    ts = rolling_average(ts,10)
    ax.plot(time,ts,color=cmap(row_indx/7),label=dc,linewidth=2)
    row_indx+=1

ax.legend(frameon=False,title='log drug \nconc')
ax.set_ylabel('Normalized fluorescence',fontsize=14)
ax.set_xlabel('Time',fontsize=14)
ax.set_title('With AB-CTX interaction calibration')
# %%
fig,ax = plt.subplots()

row_indx = 0
col = '11'

max_val = np.max(p.data['A'+col])

for row in row_list[0:-1]:
    key = row+col
    dc = dc_log[row_indx]
    norm_coeff = get_polyfit_result(dc,res_f_v_dc)
    ts = np.array(p.data[key])
    # ts = ts/norm_coeff
    ts = ts/max_val
    ts = rolling_average(ts,10)
    ax.plot(time,ts,color=cmap(row_indx/7),label=dc,linewidth=2)
    row_indx+=1

ax.legend(frameon=False,title='log drug \nconc')
ax.set_ylabel('Normalized fluorescence',fontsize=14)
ax.set_xlabel('Time',fontsize=14)
ax.set_title('Without AB-CTX interaction calibration');
# %%

od_data = p.parse_data_file(data_file,data_start='OD600')

fig,ax = plt.subplots()

row_indx = 0
col = '11'

for row in row_list[0:-1]:
    key = row+col
    dc = dc_log[row_indx]

    ts = np.array(od_data[key])

    # ts = rolling_average(ts,10)
    ax.plot(time,ts,color=cmap(row_indx/7),label=dc,linewidth=2)
    row_indx+=1

y_t = [0,1.1]
x_t = [15500,15500]
ax.plot(x_t,y_t,'--',linewidth=3)

ax.set_xlim(0,20000)
ax.set_ylim(0,0.6)
# %%

cell_count_indx = np.argwhere(np.array(time)>15500)[0][0]
od_list = []

for row in row_list:
    key = row + col
    ts = np.array(od_data[key])
    od_list.append(ts[cell_count_indx])

od_list = np.array(od_list) - np.min(od_list)

cell_count_list = [np.log10(get_cells_from_od(od)) for od in od_list]

# %% try to fit reaction kinetics for the first 3 timeseries data

exp_start = 145

col = '11'

max_val = np.max(p.data['A'+col])
row_indx = 2

ydata = []

time_t = time[exp_start:]
time_t = time_t - np.min(time_t)

fig,ax = plt.subplots()

for row in ['C']:
    key = row+col
    dc = dc_log[row_indx]
    norm_coeff = get_polyfit_result(dc,res_f_v_dc)
    ts = np.array(p.data[key])
    ts = ts/max_val
    ts = ts/norm_coeff
    ts = rolling_average(ts,10)
    ydata.append(ts[exp_start:])

    row_indx += 1


x0 = [-8,-5]
res = scipy.optimize.minimize(global_loss_known_N,x0,
                              args=(ydata,time_t,cell_count_list[0:3]),
                              method='CG',tol=10**-10,
                              options={'disp':True,'adaptive':True,'maxiter':1000})

ka = res.x[0]
kd = res.x[1]
ka = -9.5
# # kd = -5.5
row_indx=0

for yd in ydata:

    ax.plot(time_t,yd,color='black')
    y_t = kinetic_sol(time_t,[1,0],ka,kd,cell_count_list[row_indx])[:,1]
    ax.plot(time_t,y_t,color='red')
    row_indx+=1
# %%
