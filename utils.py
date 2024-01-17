import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import re
import scipy.optimize as sciopt
from scipy.integrate import odeint

def est_linear_slope(counts,
                    # dt=1, # per hour
                    dt=4,
                    time=None,
                    window=10,
                    step=5,
                    thresh=0.65,
                    debug=True,
                    title=None,
                    exclude=0,
                    return_fit=False):
    """A non-parametric method for estimating the slope of growth curves

    Args:
        OD (array-like): Raw optical density data
        window (float,optional): Window for rolling linear regression. Defaults to 10.
        step (float, optional): Step size for rolling lienar regression. Defaults to 5.
        thresh (float, optional): Percentage of max slope to include in final 
        linear regression. Defaults to 0.65.
        debug (bool, optional): If True, plots linear fit over data. Defaults to False.

    Returns:
        _type_: _description_
    """

    n_measurements = len(counts)
    if type(counts) is list:
        counts = np.array(counts)

    if time is None:
        time = np.arange(len(counts))*dt

    # remove zero values
    time = np.delete(time,np.argwhere(counts<=0))
    counts = np.delete(counts,np.argwhere(counts<=0))

    if len(counts) < n_measurements/2:
        return np.nan
    
    lnCount = np.log(counts)

    # check if the population is decreasing
    if np.mean(lnCount[-3:]) < np.mean(lnCount[:3]):
        est_negative_slope = True
    else:
        est_negative_slope = False

    slopes = []
    x_pos = []

    # calculate piecewise slopes
    for indx in range(exclude,len(lnCount)-window+step,step):
        time_t = time[indx:indx+window]
        subset = lnCount[indx:indx+window]
        fit = scipy.stats.linregress(time_t,subset)
        slopes.append(fit.slope)
        x_pos.append(indx)
    
    if est_negative_slope:
        lb = thresh*np.nanmin(slopes)
    else:
        lb = thresh*np.nanmax(slopes)

    if est_negative_slope:
        use_indices = np.argwhere(slopes<=lb)[:,0]
    else:
        use_indices = np.argwhere(slopes>=lb)[:,0]

    # if len(use_indices) == 1:
    #     if est_negative_slope:
    #         print('Warning: only one slope found')
    #         return np.nanmin(slopes)
    #     else:
    #         print('Warning: only one slope found')
    #         return np.nanmax(slopes)

    if len(use_indices) > 1:
        lin_range = x_pos[np.min(use_indices):np.max(use_indices)+window]
    else:
        lin_range = x_pos[use_indices[0]:use_indices[0]+window]

    # compute slope and plot result
    time_t = time[lin_range]
    lin_seg = lnCount[lin_range]
    fit = scipy.stats.linregress(time_t,lin_seg)
    slope = fit.slope

    count_fit = time_t*fit.slope + fit.intercept
    
    if return_fit:
        return (fit.slope,fit.intercept)

    if np.isnan(slope):
        raise Warning('Slope is NaN, adjust parameters')
        return slope

    # plot the linear regression
    if debug:
        fig,ax = plt.subplots(figsize=(4,3))
        # ax = ax_list[0]
        ax.plot(time,lnCount,linewidth=2,label='lnCount')
        ax.plot(time_t,count_fit,linewidth=2,label='Linear fit')
        # ax.plot(time[x_pos],(10000*np.array(slopes))+np.min(lnCount),label='slope')
        # ax.plot(cutoff_time,cutoff,label='Threshold')
        ax.legend(frameon=False)

        # title = title + ' ' + 'err = ' + str(round(100000*fit.stderr,2))

        ax.set_title(title)

    return slope # per hour

def plot_plate(data,plot_fit=False,
                label_well=False,row_list=None,col_list=None,dt=4,
                slope_est_options={},col_key_color=None,col_key_genotype=None,
                color_dict=None,title=None):

    if row_list is None:
        row_list = ['B','C','D','E','F','G']
    if col_list is None:
        col_list = np.arange(9) + 2
        col_list = [str(c) for c in col_list]

    ncols = len(col_list)
    nrows = len(row_list)

    fig,ax_list = plt.subplots(nrows=nrows,ncols=ncols,figsize=(12,8),
                                sharex=True,sharey=True)
    
    key0 = list(data.keys())[0]
    time = np.arange(len(data[key0]['green']))*dt

    row_indx = 0
    for row in row_list:
        col_indx = 0
        
        for col in col_list:
            ax = ax_list[row_indx,col_indx]
            key = row+col
            if col_key_color is None:
                green_data = data[key]['green']
                red_data = data[key]['red']
                ax.plot(time,np.log(green_data),color='tab:cyan',linewidth=2.5)
                ax.plot(time,np.log(red_data),color='tab:brown',linewidth=2.5)
            
            else:
                color = col_key_color[col]
                genotype = col_key_genotype[col]
                plot_color = color_dict[genotype]

                ts = data[key][color]
                time = np.arange(len(ts))*dt
                ax.plot(time,ts,color=plot_color,linewidth=3)

            if plot_fit:

                # if np.max(ts) > 10**3:
                fit = est_linear_slope(ts,time=time,return_fit=True,**slope_est_options,debug=False)
                if type(fit) == tuple:
                    fit_t = time*fit[0] + fit[1]
                    fit_t = np.e**fit_t
                    ax.plot(time,fit_t,'--',color='black',alpha=0.7)

                    ax.annotate(str(round(fit[0]*10,2)),xy=(0.05,0.8),xycoords='axes fraction',fontsize=10)

            if label_well:
                ax.set_title(key)

            ax.set_yscale('symlog',linthresh=10)

            col_indx += 1
        row_indx += 1

    if title is not None:
        fig.suptitle(title)

    fig.tight_layout()

    return fig,ax

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

def hill_fn(conc,gmax, gmin, hc, ic_50):
    
    y = gmax + ((gmin - gmax) * conc**hc) / (ic_50**hc + conc**hc)
    return y

def est_dr_params(dc,gr,sigma=None):
    gmax_est = np.max(gr)
    gmin_est = np.min(gr)
    hc_est = 1
    ic_50_est = np.median(dc)
    p0 = [gmax_est,gmin_est,hc_est,ic_50_est]
    popt,pcov = sciopt.curve_fit(hill_fn,dc,gr,p0=p0,
                            maxfev=10000,sigma=sigma)
    return popt,pcov

def cell_count_vs_time(t,N0,r_g,r_d,alpha):
    """Cell count vs time

    Args:
        t (array-like): Time
        N0 (float): Initial cell count (log)
        r_g (float): Growth rate
        r_d (float): Death rate
        alpha (float): Antibiotic lag time
    """

    N = N0 + (1/2.303)*((r_g - r_d)*t + (r_d/alpha)*(np.exp(-alpha*t)))

    return N

def growth_diffeq(N,t,K,Kss,alpha,cc):

    dydt = (K-Kss*(1-np.exp(-alpha*t)))*N*(1-N/cc)

    return dydt

def growth_sol(t,y0,K,Kss,alpha,cc):
    y = odeint(growth_diffeq,y0,t,args=(K,Kss,alpha,cc))
    return y[:,0]