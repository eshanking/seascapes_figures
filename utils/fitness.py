import numpy as np
from numpy.ma.core import get_data
import scipy
from seascapes_figures.utils import dir_manager
import matplotlib.pyplot as plt

def gen_fitness_curves(pop,conc=None):
    
    if conc is None:
        conc = np.logspace(-3,5,num=1000)
    
    n_genotype = pop.n_genotype

    fc = {}
    for g in range(n_genotype):
        f = np.zeros(len(conc))
        i = 0
        for c in conc:
            f[i] = gen_fitness(pop,g,c) - pop.death_rate
            i+=1
        fc[g] = f

    return fc

# compute fitness given a drug concentration
def gen_fitness(pop,genotype,conc,drugless_rate=None,ic50=None):        

    if drugless_rate is None:
        drugless_rate = pop.drugless_rates
    if ic50 is None:
        ic50 = pop.ic50

    # logistic equation from Ogbunugafor 2016
    conc = conc/10**6 # concentration in uM, convert to M
    c = -.6824968 # empirical curve fit
    log_eqn = lambda d,i: d/(1+np.exp((i-np.log10(conc))/c))
    if conc <= 0:
        fitness = drugless_rate[genotype]
    else:
        fitness = log_eqn(drugless_rate[genotype],ic50[genotype])

    return fitness

def logistic_equation(conc,drugless_rate,ic50):
    """
    Logistic equation from ogbunugafor et al, PLOS CB, 2016

    Parameters
    ----------
    dugless_rate : float
        Drugless growth rate of genotype.
    ic50 : float
        ic50 of genotype.
    conc : float
        Drug concentration (in Molarity (M)).
    c : float, optional
        Logistic curve steepness parameter. The default is -0.6824968.

    Returns
    -------
    f : float
        Replication rate.

    """
    c=-0.6824968
    conc = conc/10**6
    f = drugless_rate/(1+np.exp((ic50-np.log10(conc))/c))
    
    return f

def gen_static_landscape(pop,conc):

    # get final landscape and seascape
    landscape = np.zeros(pop.n_genotype)
    for kk in range(pop.n_genotype):
        landscape[kk] = gen_fitness(pop,
                                    kk,
                                    pop.static_topo_dose,
                                    pop.drugless_rates,
                                    pop.ic50)
    
    if min(landscape) == 0:
        zero_indx_land = np.argwhere(landscape==0)
        landscape_t = np.delete(landscape,zero_indx_land)
        min_landscape = min(landscape_t)
    else:
        min_landscape = min(landscape)
        zero_indx_land = []
    
    seascape = np.zeros(pop.n_genotype)
    for gen in range(pop.n_genotype):
        seascape[gen] = gen_fitness(pop,gen,conc,pop.drugless_rates,pop.ic50)
        
    if min(seascape) == 0:
        zero_indx_sea = np.argwhere(seascape==0)
        seascape_t = np.delete(seascape,zero_indx_sea)
        min_seascape = min(seascape_t)
    else:
        min_seascape = min(seascape)
        
    landscape = landscape - min_landscape
    landscape = landscape/max(landscape)
    
    rng = max(seascape) - min_seascape
    
    landscape = landscape*rng + min_seascape
    
    landscape[zero_indx_land] = 0
    return landscape

def gen_digital_seascape(pop,conc,gen,min_fitness=0):
    if pop.mic_estimate is not None:
        mic = est_mic(pop,gen,Kmic=pop.mic_estimate)
    else:
        mic = est_mic(pop,gen,growth_rate=pop.death_rate)
    
    if conc >= mic:
        fitness = min_fitness
    else:
        fitness = pop.drugless_rates[gen]
    return fitness

def gen_fit_land(pop,conc,mode=None):
    
    fit_land = np.zeros(pop.n_genotype)
            
    if pop.fitness_data == 'manual' or mode=='manual':
        fit_land = pop.landscape_data/pop.doubling_time
   
    else:
        
        if pop.static_topology:
            fit_land = gen_static_landscape(pop,conc)
            
        if pop.digital_seascape:
            for kk in range(pop.n_genotype):
                fit_land[kk] = gen_digital_seascape(pop, conc, kk)
            
        else:
            for kk in range(pop.n_genotype):
                fit_land[kk] = gen_fitness(pop,
                                           kk,
                                           conc,
                                           pop.drugless_rates,
                                           pop.ic50)/pop.doubling_time
    
    return fit_land

# Generate fitness landscape for use in the abm method
# Private to avoid confusion with gen_fit_land
def gen_fl_for_abm(pop,conc,counts):
    
    fit_land = gen_fit_land(pop,conc)
    
    # # takes the landscape at the max dose and scales the replication rate
    # # according to drug concentration
    # if pop.static_landscape:
    #     # max_fitness = max(fit_land)
    #     # fit_land = pop.gen_fit_land(pop.max_dose)
    #     # fit_land = fit_land*max_fitness/max(fit_land)
    #     fit_land = gen_fit_land(pop,conc)
    
    # if pop.static_topology:
    #     fit_land = gen_fit_land(pop,conc)
    
    # Scale division rates based on carrying capacity
    if pop.carrying_cap:
        division_scale = 1-np.sum(counts)/pop.max_cells
        if counts.sum()>pop.max_cells:
            division_scale = 0
    else:
        division_scale = 1
    
    fit_land = fit_land*division_scale
    
    return fit_land

def gen_random_seascape(n_allele,
                        drugless_limits=[1,1.5],
                        ic50_limits=[-6.5,-1.5]):
    
    n_genotype = 2**n_allele
    
    drugless_rates = np.random.uniform(min(drugless_limits),
                                      max(drugless_limits),
                                      n_genotype)
    
    ic50 = np.random.uniform(min(ic50_limits),
                             max(ic50_limits),
                             n_genotype)
    
    return drugless_rates,ic50

def randomize_seascape(pop,
                       drugless_limits=[1,1.5],
                       ic50_limits=[-6.5,-1.5]):
    
    n_genotype = pop.n_genotype
    
    pop.drugless_rates = np.random.uniform(min(drugless_limits),
                                           max(drugless_limits),
                                           n_genotype)

    pop.ic50 = np.random.uniform(min(ic50_limits),
                                 max(ic50_limits),
                                 n_genotype)
    
def fit_logistic_curve(xdata,ydata):
    from scipy.optimize import curve_fit
    
    popt,var = curve_fit(logistic_equation,xdata,ydata)
    
    return popt

def gen_null_seascape(pop,conc):

    
    landscape = gen_fit_land(pop,conc)
    start_rates = gen_fit_land(pop,10**-3)
    final_rates = gen_fit_land(pop,10**5)
    # mid_rates = gen_fit_land(pop,10**1)
    
    start_points = scale_and_ignore_zeros(landscape,start_rates)
    end_points = scale_and_ignore_zeros(landscape,final_rates)
    # mid_points = scale_and_ignore_zeros(landscape,mid_rates)
    mid_points = landscape
    
    xdata = [10**-3,conc,10**5]
    
    ic50_new = []
    drugless_rates_new = []
    
    for genotype in range(len(landscape)):
        ydata = [start_points[genotype],
                 mid_points[genotype],
                 end_points[genotype]]
        params = fit_logistic_curve(xdata,ydata)
        ic50_new.append(params[1])
        drugless_rates_new.append(params[0])
    # find the null landscape drugless rates
    
    drugless_rates_new = scale_and_ignore_zeros(drugless_rates_new,
                                                pop.drugless_rates)
    
    return drugless_rates_new,ic50_new

def scale_and_ignore_zeros(data,target):
    """
    Scale data to range of target while ignoring the zero values in data and
    target.

    Parameters
    ----------
    data : numpy array
        Data to be scaled to the range of target.
    target : numpy array
        Target data range.

    Returns
    -------
    scaled_data : numpy array
        Scaled data to range of target. Zero values in data are set to zero
        in scaled_data and zero values in target are not used to calculate
        range.

    """
    # make sure inputs are numpy arrays
    
    if not isinstance(data,np.ndarray):
        data=np.array(data)
    if not isinstance(target,np.ndarray):
        target=np.array(target)
    
    if min(data) == 0:
        zero_indx_data = np.argwhere(data==0)
        data_t = np.delete(data,zero_indx_data)
        min_data = min(data_t)
    else:
        min_data = min(data)
        zero_indx_data = []
        
    if min(target) == 0:
        zero_indx_target = np.argwhere(target==0)
        target_t = np.delete(target,zero_indx_target)
        min_target = min(target_t)
    else:
        min_target = min(target)
        
    data = data - min_data
    data = data/max(data)

    rng = max(target) - min_target
    
    scaled_data = data*rng + min_target

    scaled_data[zero_indx_data] = 0
    
    return scaled_data

def est_mic(pop,gen,Kmic=None,growth_rate=None):
    """
    est_mic: estimates the mic based on a given Kmic (ratio of growth rate to 
    max growth rate at MIC) or based on a given growth rate.

    Parameters
    ----------
    pop : population class object
        
    gen : int
        Genotype under consideration.
    Kmic : float, optional
        Ratio of growth rate to max growth rate at MIC. The default is None.
    growth_rate : float, optional
        Growth rate at MIC. The default is None.

    Raises
    ------
    Exception
        Function requires Kmic OR growth_rate to calculate MIC.

    Returns
    -------
    mic : float
        MIC at a given growth rate or Kmic.

    """
    
    if Kmic is None:
        if growth_rate is None:
            raise Exception('Need a growth rate or Kmic threshold to estimate mic.')
        else:
            Kmic = growth_rate/pop.drugless_rates[gen]
    c=-0.6824968
    mic = 10**(pop.ic50[gen]+6 - c*np.log((1/Kmic)-1))
    return mic

def get_background_keys(df):

    # row A, row H, col 1, and col 12

    k = df.keys()

    k = k[2:]
    bg_keys = [y for y in k if int(y[1:]) == 1] # col 1
    bg_keys = bg_keys + [y for y in k if (int(y[1:]) == 12 and y not in bg_keys)]
    bg_keys = bg_keys + [y for y in k if (y[0] == 'A' and y not in bg_keys)]
    bg_keys = bg_keys + [y for y in k if (y[0] == 'H' and y not in bg_keys)]

    return bg_keys

def get_data_keys(df):

    bg_keys = get_background_keys(df)

    data_keys = [k for k in df.keys() if k not in bg_keys]
    data_keys = data_keys[2:]

    return data_keys

def estimate_background(df):

    bg_keys = get_background_keys(df)
    s = 0

    for key in bg_keys:
        s += np.average(df[key])

    bg = s/len(bg_keys)

    return bg


def subtract_background(df):

    bg = estimate_background(df)
    datakeys = df.keys()[2:]

    if (df.values < bg).any():
        bg = np.min(df.values)
        

    for key in datakeys:
        df[key] = df[key] - bg

    return df

def get_growth_rate_data(data_path):
    
    df = dir_manager.load_growth_rate_data(data_path)
    df = subtract_background(df)

    return df

def get_growth_rates_from_df(df,carrying_cap=None):
    
    """Extracts the growth rates from timeseries growth data in a dataframe

    Arguments:
        df: pandas DataFrame
            data frame containing growth rate data

    Returns:
        growth_rates: dict
            dictionary of growth rates for each experimental condition
    """

    growth_rates = {}

    # assumes outer row is background data
    data_keys = get_data_keys(df)
    time = df['Time [s]']

    for k in data_keys:
        growth_rates[k] = est_growth_rate(df[k],t=time,carrying_cap=carrying_cap)
        # cur_genotype = k[0]
        # if cur_genotype == prev_genotype:
        #     growth_rates[k] = est_growth_rate(df[k],t=time,carrying_cap=cc)
        # else:
        #     # get new carrying capacity
        #     genotype_keys = [k for k in data_keys if k[0]==prev_genotype]


    return growth_rates

def gen_seascape_library(pop):

    seascape_lib = {}
    seascape_lib['drug_conc'] = pop.seascape_drug_conc # add drug concentration data

    library = ['B','C','D','E','F'] # library of column names
    genotype = 0

    for df in pop.growth_rate_data: # for each plate
        
        growth_rates = get_growth_rates_from_df(df,carrying_cap=pop.max_od) # estimate the growth rate for each set of timeseries data

        for l in library:
            gr_vect = []
            i = 2
            for c in pop.seascape_drug_conc:
                key = l + str(i) # each key is a seperate column in the plate reader data,looping over i holds genotype steady and increments drug concentration
                gr_vect.append(growth_rates[key])
                i+=1
            # print(len(gr_vect))
            seascape_lib[str(genotype)] = gr_vect
            genotype += 1
    
    # df3 = pop.growth_rate_data[2]
    # print(genotype)
    l = 'G'
    i=2
    gr_vect = []
    for c in pop.seascape_drug_conc:
        key = l + str(i)
        # print(key)
        gr_vect.append(growth_rates[key])
        i+=1
    
    seascape_lib[str(genotype)] = gr_vect

    return seascape_lib

def est_growth_rate(growth_curve,t=None,debug=False,carrying_cap=None):

    """Estimates microbial growth rate from OD growth rate curves

    Arguments:
        growth_rate: list of floats or numpy array
            Growth rate data
        t (optional): list of floats or numpy array. 
            Time vector corresponding to growth rate data. If t is not given, algorithm 
            assumes each time step is 1 second.
        degug (optional): boolean
            Default False. If True, displays plots useful for debugging

    Returns:
        float: Growth rate in units of 1/s
    """

    if t is None:
        t = np.arange(len(growth_curve))

    if not isinstance(growth_curve,np.ndarray):
        growth_curve = np.array(growth_curve)
    
    # normalize to range (0,1)
    growth_curve = growth_curve - min(growth_curve)
    if carrying_cap is not None:
        growth_curve = growth_curve/carrying_cap
    else:
        growth_curve = growth_curve/max(growth_curve)

    # method 1
    # compute derivative
    dpdt = np.zeros(len(growth_curve)-1)
    r = np.zeros(len(dpdt))
    for i in range(len(growth_curve)-1):
        p = growth_curve[i]
        dp = growth_curve[i+1] - growth_curve[i]
        dt = t[i+1] - t[i]
        dpdt[i] = dp/dt
        # compute r
        if p == 0 or 1-p == 0:
            r[i] = 0
        else:
            r[i] = dpdt[i]/(1/(p*(1-p)))
    
    r1 = np.max(r)
    if debug:
        fig,ax = plt.subplots(nrows=3)
        ax[0].hist(r)
        ax[1].plot(t,growth_curve)
        ax[2].plot(t[:-1],r)
        ax[0].set_title('Growth rate histogram')
        ax[0].set_xlabel('$r$')
        ax[1].set_title('Normalized growth curve')
        ax[1].set_xlabel('t (s)')
        ax[2].set_title('Growth rate over time')
        ax[2].set_xlabel('t (s)')
        plt.tight_layout()
        
    return r1

def get_max_od(pop):

    max_od = 0
    for df in pop.growth_rate_data:
        data_keys = get_data_keys(df)
        for key in data_keys:
            if max_od < max(df[key]):
                max_od = max(df[key])

    return max_od