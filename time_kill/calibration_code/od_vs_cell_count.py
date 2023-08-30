from fears.utils import AutoRate
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline
import pickle

od_data = 'calibration_05162023/EK_single_OD600_20230516_103116_with_lid.xlsx'
p_od = AutoRate.Plate(od_data,mode='single_measurement')
od_data = p_od.od_data_to_dict(p_od.data)

row_list = ['A','B','C','D','E','F','G','H']
bg_cols = [1,12]
bg_cols = [str(c) for c in bg_cols]

# estimate od background    

bg_est = 0
indx = 0
for row in row_list:
    for col in bg_cols:
        key = row + col
        bg_est += od_data[key]
        indx+=1

bg_est = bg_est/indx

# subtract background

col_list = np.arange(12) + 1
col_list = [str(c) for c in col_list]

od_t = np.zeros((len(row_list)-2,len(col_list)-2))
row_indx = 0
for row in row_list[1:-1]:
    col_indx = 0
    for col in col_list[1:-1]:
        key = row+col
        od_data[key] = od_data[key] - bg_est
        od_t[row_indx,col_indx] = od_data[key]
        col_indx += 1
    row_indx += 1

od_avg = np.mean(od_t,axis=0)
od_std = np.std(od_t,axis=0)
od_err = od_std/np.sqrt(len(row_list)-2)


# plot OD vs dilution
dilution_str = []
dilution = []

for i in range(len(col_list)-2):
    dilution_str.append('1/'+str(i+1))
    dilution.append(1/(i+1))

fig,ax = plt.subplots()

K = 9*10**5 # cells/uL
cell_count = K*np.array(dilution)

ax.errorbar(od_avg,cell_count,xerr=od_err,fmt='o',color='k',capsize=5,markersize=7)
# ax.set_yscale('log')
ax.set_xscale('log')
ax.set_yscale('log')

def logistic_fn(x,a,b,c):
    return a/(1+np.exp(-b*(x-c)))

p0 = [np.max(cell_count),20,np.mean(od_avg)]
bounds = [[0,0,0],[10,10,10]]
popt,pcov = curve_fit(logistic_fn,od_avg,cell_count,p0=p0)

xfit = np.linspace(np.min(od_avg),np.max(od_avg),100)
yfit = logistic_fn(xfit,*popt)
# ax.plot(xfit,yfit,'k--',alpha=0.5)

# try spline fit
# reverse data
od_avg = od_avg[::-1]
cell_count = cell_count[::-1]
spl = UnivariateSpline(od_avg,cell_count,s=12)
yfit = spl(xfit)
ax.plot(xfit,yfit,'r--')

# pickle the spline for later use

with open('od_cell_count_with_lid_spline.pkl', 'wb') as file:
    pickle.dump(spl, file)

# test the pickled spline

with open('od_cell_count_with_lid_spline.pkl', 'rb') as file:
    spl = pickle.load(file)

yfit = spl(xfit)
ax.plot(xfit,yfit,'b--')