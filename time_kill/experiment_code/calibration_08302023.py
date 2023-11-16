#%%
from fears.utils import AutoRate
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import stats
import scipy.interpolate as interp
import scipy.optimize as sciopt

def run():
    ab_plate_path = '../calibration_data/calibration_08302023/EK_single_AB_constant_gain_20230830_122655.xlsx'
    od_plate_path = '../calibration_data/calibration_08302023/EK_single_OD600_20230830_114542.xlsx'
    p_ab = AutoRate.Plate(ab_plate_path,mode='single_measurement')
    p_od = AutoRate.Plate(od_plate_path,mode='single_measurement')
    od_data = p_od.od_data_to_dict(p_od.data)
    row_list = ['A','B','C','D','E','F','G','H']
    bg_cols = [1,11,12]
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

#%% Estimate cell count from OD
    row_indx = 0

    # refactor data into 2d array
    od_data_t = np.zeros((len(row_list)-2,len(col_list)-2))
    for row in row_list[1:-1]:
        col_indx = 0
        for col in col_list[1:-1]:
            key = row + col
            od_data_t[row_indx,col_indx] = od_data[key]
            col_indx += 1
        row_indx += 1

    # calculate mean and standard error
    od_avg = np.mean(od_data_t,axis=0)
    od_std = np.std(od_data_t,axis=0)
    od_err = od_std/np.sqrt(len(row_list)-2)

    # fig,ax = plt.subplots()

    # generate dilutions based on 2x dilution scheme

    dilutions_str = []
    dilutions = []
    for i in range(len(col_list)-2):
        dilutions_str.append(str(2**i) + 'x')
        dilutions.append(1/(2**i))

    dilutions[-1] = 0
    dilutions_str[-1] = '0x'


    def od_to_cells(od):
        """Background subtracted od to cell count

        Args:
            od (_type_): _description_

        Returns:
            _type_: _description_
        """
        m = 1
        b = 13.52
        res = []
        for o in od:
            if o < 0:
                res.append(0)
            else:
                res.append(np.exp(m*np.log(o) + b))
        return np.array(res)

    cell_count = od_to_cells(od_avg) # estimated from od data
    cell_count_err = od_to_cells(od_err)

    cell_count_est = [d*cell_count[0] for d in dilutions] # estimated from dilution scheme
    cell_count_est[-1] = 0

    # get the fluorescence data

    fluor_data = p_ab.od_data_to_dict(p_ab.data)

    cmap = mpl.colormaps['viridis']

    # estimate background
    bg_est = 0
    indx = 0
    for row in row_list:
        col = '11' # no cells in this column
        key = row + col
        bg_est += fluor_data[key]
        indx+=1
    bg_est = bg_est/indx

    fluor_data_t = np.zeros((len(row_list)-2,len(col_list)-2))
    fluor_data_std = np.zeros((len(row_list)-2,len(col_list)-2))
    row_indx = 0
    for row in row_list[1:-1]:
        col_indx = 0
        for col in col_list[1:-1]:
            key = row + col
            fluor_data_t[row_indx,col_indx] = fluor_data[key] - bg_est
            col_indx += 1
        row_indx += 1

    fluor_avg = np.mean(fluor_data_t,axis=0) # background substrated fluorescence
    fluor_err = np.std(fluor_data_t,axis=0)/np.sqrt(len(row_list)-2)

    cell_count_norm = cell_count[1:]/cell_count[1]
    fluor_norm = fluor_avg[1:]/fluor_avg[1]

    s = 0
    k = 2

    # sort data
    sort_indx = np.argsort(fluor_norm)
    fluor_norm = fluor_norm[sort_indx]
    cell_count_norm = cell_count_norm[sort_indx]

    spline = interp.UnivariateSpline(fluor_norm,cell_count_norm,k=k,s=s)

    def rfu30_to_cell_count(rfu30,background_subtracted=True):
        """Convert RFU30 to cell count

        Args:
            rfu30 (float): RFU30 value
            background_subtracted (bool, optional): If the RFU30 value is background subtracted. Defaults to True.

        Returns:
            float: cell count
        """
        if background_subtracted == False:
            rfu30 = rfu30 - bg_est
        rfu30_norm = rfu30/fluor_avg[1]
        return spline(rfu30_norm)*cell_count[1]
    
    return rfu30_to_cell_count


# def run():
#     # ab_plate_path = 'calibration_05172023/EK_single_AB_constant_gain_20230517_161938_30_min.xlsx'
#     # od_plate_path = 'calibration_05172023/EK_single_OD600_20230517_154009_no_lid.xlsx'
#     ab_plate_path = '../calibration_data/calibration_08302023/EK_single_AB_constant_gain_20230830_122655.xlsx'
#     od_plate_path = '../calibration_data/calibration_08302023/EK_single_OD600_20230830_114542.xlsx'
#     p_ab = AutoRate.Plate(ab_plate_path,mode='single_measurement')
#     p_od = AutoRate.Plate(od_plate_path,mode='single_measurement')
#     od_data = p_od.od_data_to_dict(p_od.data)

#     #%% estimate od background

#     row_list = ['A','B','C','D','E','F','G','H']
#     bg_cols = [1,11,12]
#     bg_cols = [str(c) for c in bg_cols]

#     bg_est = 0
#     indx = 0
#     for row in row_list:
#         for col in bg_cols:
#             key = row + col
#             bg_est += od_data[key]
#             indx+=1

#     bg_est = bg_est/indx

#     col_list = np.arange(12) + 1
#     col_list = [str(c) for c in col_list]

#     for row in row_list:
#         for col in col_list:
#             key = row+col
#             od_data[key] = od_data[key] - bg_est

#     #%% Estimate cell count from OD
#     row_indx = 0

#     od_data_t = np.zeros((len(row_list)-2,len(col_list)-2))
#     for row in row_list[1:-1]:
#         col_indx = 0
#         for col in col_list[1:-1]:
#             key = row + col
#             od_data_t[row_indx,col_indx] = od_data[key]
#             col_indx += 1
#         row_indx += 1

#     od_avg = np.mean(od_data_t,axis=0)
#     od_std = np.std(od_data_t,axis=0)
#     od_err = od_std/np.sqrt(len(row_list)-2)

#     # fig,ax = plt.subplots()

#     dilutions_str = []
#     dilutions = []
#     for i in range(len(col_list)-2):
#         dilutions_str.append(str(2**i) + 'x')
#         dilutions.append(1/(2**i))

#     dilutions[-1] = 0
#     dilutions_str[-1] = '0x'

#     # ax.errorbar(dilutions_str[1:],od_avg[1:],yerr=od_err[1:],fmt='x',capsize=5,color='k')
#     # ax.errorbar(dilutions_str,od_avg,yerr=od_err,fmt='x',capsize=5,color='k')
#     # ax.set_yscale('log')

#     #%%

#     def od_to_cells(od):
#         """Background subtracted od to cell count

#         Args:
#             od (_type_): _description_

#         Returns:
#             _type_: _description_
#         """
#         m = 1
#         b = 13.52
#         res = []
#         for o in od:
#             if o < 0:
#                 res.append(0)
#             else:
#                 res.append(np.exp(m*np.log(o) + b))
#         return np.array(res)

#     cell_count = od_to_cells(od_avg)
#     cell_count_err = od_to_cells(od_err)

#     cell_count_est = [d*cell_count[0] for d in dilutions]
#     cell_count_est[-1] = 0

#     # fig,ax = plt.subplots()

#     # ax.plot(cell_count_est,cell_count,'o',color='k')
#     # ax.set_xscale('log')
#     # ax.set_yscale('log')

#     # (ymin,ymax) = ax.get_ylim()
#     # (xmin,xmax) = ax.get_xlim()

#     # ax.set_ylim([np.min([xmin,ymin]),np.max([xmax,ymax])])  
#     # ax.set_xlim([np.min([xmin,ymin]),np.max([xmax,ymax])]) 

#     # ax.set_xlabel('Dilution estimation method',fontsize=14)
#     # ax.set_ylabel('OD600 estimation method',fontsize=14)

#     # # linear fit
#     # res = stats.linregress(cell_count_est,cell_count)

#     # xfit = np.linspace(np.min(cell_count_est),np.max(cell_count_est),100)
#     # yfit = res.slope*xfit + res.intercept
#     # ax.plot(xfit,yfit,'--',color='k',label='linear fit')


#     # %%
#     # fig,ax_list = plt.subplots(nrows=4,figsize=(4,8))
#     fluor_data = p_ab.od_data_to_dict(p_ab.data)

#     cmap = mpl.colormaps['viridis']

#     fluor_data_t = np.zeros((len(row_list)-2,len(col_list)-2))
#     fluor_data_std = np.zeros((len(row_list)-2,len(col_list)-2))
#     row_indx = 0
#     for row in row_list[1:-1]:
#         col_indx = 0
#         for col in col_list[1:-1]:
#             key = row + col
#             fluor_data_t[row_indx,col_indx] = fluor_data[key]
#             col_indx += 1
#         row_indx += 1

#     fluor_avg = np.mean(fluor_data_t,axis=0)
#     fluor_err = np.std(fluor_data_t,axis=0)/np.sqrt(len(row_list)-2)

#     # fig,ax = plt.subplots()

#     # ax.errorbar(fluor_avg[1:],cell_count[1:],xerr=fluor_err[1:],fmt='o',capsize=5,color='k')

#     # # ax.set_xscale('log')
#     # ax.set_yscale('log')
#     # ax.spines['right'].set_visible(False)
#     # ax.spines['top'].set_visible(False)

#     # ax.set_xlabel('Fluorescence (RFU$_{30}$)',fontsize=12)
#     # ax.set_ylabel('Cell Count (cells/$\mu$L)',fontsize=12)
#     #%%
#     # fig,ax = plt.subplots()

#     # ax.errorbar(fluor_avg[1:],cell_count[1:],xerr=fluor_err[1:],fmt='o',capsize=5,color='k')

#     # fluor_err_norm = fluor_err[1:]/np.max(fluor_avg[1:])
#     fluor_norm = fluor_data_t/np.max(fluor_data_t)
#     fluor_norm_avg = np.mean(fluor_norm,axis=0)
#     fluor_norm_err = np.std(fluor_norm,axis=0)/np.sqrt(len(row_list)-2)

#     cell_count_norm = cell_count/np.max(cell_count)

#     # remove first point
#     fluor_norm_avg = fluor_norm_avg[1:]
#     fluor_norm_err = fluor_norm_err[1:]

#     cell_count_norm = cell_count_norm[1:]

#     # def inv_logistic(y,y0,k,x50,l):
#     #     return x50/((y0/(y-l) - 1)**k)

#     # def inv_logistic(x,xmax,k,y50,l):
#     #     return y50/(xmax/(x-l) - 1)**k

#     # def power_law_func(X, a, b, c, d):
#     #     # X = X - l
#     #     return a * np.power(X, b) + c * np.power(X, d)

#     # p0 = [1.1,1,0.1,0.1]
#     # popt_unweighted,pcov = sciopt.curve_fit(power_law_func,fluor_norm_avg,cell_count_norm,
#     #                                       p0=p0,maxfev=10000)
#     # # # popt,pcov = sciopt.curve_fit(expon,RFU60_avg,cell_count[1:],p0=[1,1,1])
#     # # popt[-1] = 
#     # xfit = np.linspace(np.min(fluor_norm_avg),np.max(fluor_norm_avg),100)
#     # yfit = power_law_func(xfit,*popt_unweighted)

#     # ax.plot(xfit*np.max(fluor_data_t),yfit*np.max(cell_count),label='power law unweighted',linewidth=2)

#     # ax.set_yscale('log')

#     # # now try weighted fit

#     # popt_weighted,pcov = sciopt.curve_fit(power_law_func,fluor_norm_avg,cell_count_norm,
#     #                                       p0=p0,sigma=fluor_norm_err,maxfev=10000)

#     # xfit = np.linspace(np.min(fluor_norm_avg),np.max(fluor_norm_avg),100)
#     # yfit = power_law_func(xfit,*popt_weighted)

#     # ax.plot(xfit*np.max(fluor_data_t),yfit*np.max(cell_count),label='power law weighted',linewidth=2)
#     # ax.legend(frameon=False,fontsize=12)
#     # ax.spines['right'].set_visible(False)
#     # ax.spines['top'].set_visible(False)

#     # ax.set_xlabel('Fluorescence (RFU$_{30}$)',fontsize=12)
#     # ax.set_ylabel('Cell Count (cells/$\mu$L)',fontsize=12)
#     #%%
#     # fig,ax = plt.subplots()
#     # now try spline

#     s = 0
#     k = 2

#     # sort data
#     sort_indx = np.argsort(fluor_norm_avg)
#     fluor_norm_avg = fluor_norm_avg[sort_indx]
#     cell_count_norm = cell_count_norm[sort_indx]

#     spline = interp.UnivariateSpline(fluor_norm_avg,cell_count_norm,k=k,s=s)
#     xfit = np.linspace(np.min(fluor_norm_avg),np.max(fluor_norm_avg),100)
#     yfit = spline(xfit)

#     # ax.plot(xfit*np.max(fluor_data_t),yfit*np.max(cell_count),label='spline',linewidth=2)

#     # # ax.legend(frameon=False,fontsize=12)
#     # ax.spines['right'].set_visible(False)
#     # ax.spines['top'].set_visible(False)

#     # ax.set_xlabel('Fluorescence (RFU$_{30}$)',fontsize=12)
#     # ax.set_ylabel('Cell Count (cells/$\mu$L)',fontsize=12)
#     # ax.errorbar(fluor_avg[1:],cell_count[1:],xerr=fluor_err[1:],yerr=cell_count_err[1:],fmt='o',capsize=5,color='k')

#     # ax.set_yscale('symlog',linthresh=4000)
#     # ax.set_xscale('log')
#     # ax.legend(frameon=False,fontsize=12)


#     #%%

#     def rfu30_to_cell_count(rfu30,spline):
#         rfu_norm = np.max(fluor_data_t)
#         cc_norm = np.max(cell_count)
#         return cc_norm*spline(rfu30/rfu_norm)

#     #%% pickle objects

#     # with open('../spline.pkl','wb') as f:
#     #     pickle.dump(spline,f)

#     # with open('../rfu30_to_cell_count.pkl','wb') as f:
#     #     pickle.dump(rfu30_to_cell_count,f)
#     # %% define dynamic range

#     cell_count_range = cell_count[1:]
#     # dilution_range = dilutions_str[1:]
#     dilution_range = dilutions[1:]
#     fluor_range = fluor_avg[1:]

#     # fig,ax_list = plt.subplots(ncols=2,figsize=(8,4))

#     # ax = ax_list[0]
#     # x = np.arange(len(dilution_range))
#     # ax.plot(x,cell_count_range,'o',color='k')

#     # # ax.set_xscale('symlog',linthresh=0.0001)
#     # ax.set_yscale('symlog',linthresh=4000)

#     # ax.set_xlabel('Dilution',fontsize=12)
#     # ax.set_ylabel('Cell Count (cells/$\mu$L)',fontsize=12)

#     # ax = ax_list[1]
#     # ax.plot(x,fluor_range,'o',color='k')

#     # ax.set_yscale('log')
#     # ax.set_xlabel('Dilution',fontsize=12)
#     # ax.set_ylabel('Fluorescence (RFU$_{30}$)',fontsize=12)

#     # dilution_plot = ['2$^{1}$','2$^{2}$','2$^{3}$','2$^{4}$','2$^{5}$','2$^{6}$','2$^{7}$','2$^{8}$','0x']

#     # for ax in ax_list:
#     #     ax.spines['right'].set_visible(False)
#     #     ax.spines['top'].set_visible(False)
#     #     ax.set_xticks(x)
#     #     ax.set_xticklabels(dilution_plot,rotation=45)

#     # fig.tight_layout()

#     #%%

#     # fig,ax = plt.subplots()

#     fluor_log = np.log10(fluor_range)

#     cell_count_plot = np.array(dilution_range)*10**6

#     # ax.plot(fluor_log,cell_count_plot,'x',color='k',markersize=8,mew=2,label='data')

#     #  fit fluor vs dilution to exponential

#     def expon(x,a,r,k):
#         return a*np.exp(r*x) - k

#     p0 = [0.5,0,0]
#     popt,pcov = sciopt.curve_fit(expon,fluor_log,dilution_range,p0=p0)

#     # xfit = np.linspace(np.min(fluor_log),np.max(fluor_log),100)
#     # yfit = expon(xfit,*popt) * 10**6

#     # ax.plot(xfit,yfit,'--',color='k',label='exponential fit',linewidth=2)

#     # ax.legend(frameon=False,fontsize=12)

#     # ax.spines['right'].set_visible(False)
#     # ax.spines['top'].set_visible(False)

#     # ax.set_xlabel('Log$_{10}$ fluorescence',fontsize=12)
#     # ax.set_yscale('symlog',linthresh=3000)
#     # # ax.set_ylabel('Dilution',fontsize=12)
#     # ax.set_ylabel('Cell count (cells/$\mu$L)',fontsize=12)

#     # ax.tick_params(axis='both',which='major',labelsize=12)

#     # fig.savefig('fluor_vs_dilution.pdf',bbox_inches='tight')

#     # ax.annotate('Min fluor = {}'.format(np.min(fluor_range)),(0.05,0.4),xycoords='axes fraction')
#     # ax.annotate('Max fluor = {}'.format(np.round(np.max(fluor_range))),(0.05,0.3),xycoords='axes fraction')
#     # ax.set_yscale('log')

#     # %%
#     def rfu_to_dilution(rfu):
#         return expon(np.log10(rfu),*popt)

#     # with open('../rfu_to_dilution.pkl','wb') as f:
#     #     pickle.dump(rfu_to_dilution,f)
#     return rfu_to_dilution

# # %%
# if __name__ == "__main__":
#     rfu_to_dilution = run()