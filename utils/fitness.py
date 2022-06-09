import numpy as np
from numpy.ma.core import get_data
import scipy
from seascapes_figures.utils import dir_manager, results_manager
import matplotlib.pyplot as plt
import os

class Fitness:
    
    def __init__(self):
        return

    def gen_fitness_curves(self,pop=None,conc=None):
        
        if pop is None:
            pop = self

        if conc is None:
            conc = np.logspace(-3,5,num=1000)
        
        n_genotype = pop.n_genotype

        fc = {}
        for g in range(n_genotype):
            f = np.zeros(len(conc))
            i = 0
            for c in conc:
                f[i] = self.gen_fitness(g,c) - pop.death_rate
                i+=1
            fc[g] = f

        return fc

    # compute fitness given a drug concentration
    def gen_fitness(self,genotype,conc,drugless_rate=None,ic50=None,pop=None):        

        if pop is None:
            pop = self

        if pop.fitness_data == 'estimate':
            fitness = self.sl_to_fitness(genotype,conc)
            fitness = fitness*(60**2)

        else:
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

    def logistic_equation(self,conc,drugless_rate,ic50):
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

    def gen_static_landscape(self,conc,pop=None):
        if pop is None:
            pop = self
        # get final landscape and seascape
        landscape = np.zeros(pop.n_genotype)
        for kk in range(pop.n_genotype):
            landscape[kk] = self.gen_fitness(pop,
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
            seascape[gen] = self.gen_fitness(pop,gen,conc,pop.drugless_rates,pop.ic50)
            
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

    # def gen_digital_seascape(self,conc,gen,min_fitness=0,pop=None):

    #     if pop is None:
    #         pop = self

    #     if pop.mic_estimate is not None:
    #         mic = self.est_mic(pop,gen,Kmic=pop.mic_estimate)
    #     else:
    #         mic = self.est_mic(pop,gen,growth_rate=pop.death_rate)
        
    #     if conc >= mic:
    #         fitness = min_fitness
    #     else:
    #         fitness = pop.drugless_rates[gen]
    #     return fitness

    def gen_fit_land(self,conc,mode=None,pop=None):

        if pop is None:
            pop = self

        fit_land = np.zeros(pop.n_genotype)
                
        if pop.fitness_data == 'manual' or mode=='manual':
            fit_land = pop.landscape_data/pop.doubling_time
    
        else:
            
            if pop.static_topology:
                fit_land = self.gen_static_landscape(pop,conc)
                
            if pop.digital_seascape:
                for kk in range(pop.n_genotype):
                    fit_land[kk] = self.gen_digital_seascape(pop, conc, kk)
                
            else:
                for kk in range(pop.n_genotype):
                    fit_land[kk] = self.gen_fitness(kk,
                                            conc)/pop.doubling_time
        
        return fit_land

    # Generate fitness landscape for use in the abm method
    # Private to avoid confusion with gen_fit_land
    def gen_fl_for_abm(self,conc,counts,pop=None):

        if pop is None:
            pop = self

        fit_land = self.gen_fit_land(conc)
        
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

    def gen_random_seascape(self,n_allele,
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

    def randomize_seascape(self,pop=None,
                        drugless_limits=[1,1.5],
                        ic50_limits=[-6.5,-1.5]):
        
        if pop is None:
            pop = self

        n_genotype = pop.n_genotype
        
        pop.drugless_rates = np.random.uniform(min(drugless_limits),
                                            max(drugless_limits),
                                            n_genotype)

        pop.ic50 = np.random.uniform(min(ic50_limits),
                                    max(ic50_limits),
                                    n_genotype)
        
    def fit_logistic_curve(self,xdata,ydata):
        from scipy.optimize import curve_fit
        
        popt,var = curve_fit(self.logistic_equation,xdata,ydata)
        
        return popt

    def gen_null_seascape(self,conc,pop=None):

        if pop is None:
            pop = self

        landscape = self.gen_fit_land(conc,pop=pop)
        start_rates = self.gen_fit_land(10**-3,pop=pop)
        final_rates = self.gen_fit_land(10**5,pop=pop)
        # mid_rates = gen_fit_land(pop,10**1)
        
        start_points = self.scale_and_ignore_zeros(landscape,start_rates)
        end_points = self.scale_and_ignore_zeros(landscape,final_rates)
        # mid_points = scale_and_ignore_zeros(landscape,mid_rates)
        mid_points = landscape
        
        xdata = [10**-3,conc,10**5]
        
        ic50_new = []
        drugless_rates_new = []
        
        for genotype in range(len(landscape)):
            ydata = [start_points[genotype],
                    mid_points[genotype],
                    end_points[genotype]]
            params = self.fit_logistic_curve(xdata,ydata)
            ic50_new.append(params[1])
            drugless_rates_new.append(params[0])
        # find the null landscape drugless rates
        
        drugless_rates_new = self.scale_and_ignore_zeros(drugless_rates_new,
                                                    pop.drugless_rates)
        
        return drugless_rates_new,ic50_new

    def scale_and_ignore_zeros(self,data,target):
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

    def est_mic(self,gen,Kmic=None,growth_rate=None,pop=None):
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
        
        if pop is None:
            pop = self

        if Kmic is None:
            if growth_rate is None:
                raise Exception('Need a growth rate or Kmic threshold to estimate mic.')
            else:
                Kmic = growth_rate/pop.drugless_rates[gen]
        c=-0.6824968
        mic = 10**(pop.ic50[gen]+6 - c*np.log((1/Kmic)-1))
        return mic

################################################################
# Code for estimating fitness seascapes from OD data


    def get_background_keys(self,df):
        """Gets the dataframe keys for the background assuming a 1-well moat

        Args:
            df (pandas dataframe): dataframe containing raw OD data

        Returns:
            list: list of background keys
        """
        # row A, row H, col 1, and col 12

        k = df.keys()

        k = k[2:]
        bg_keys = [y for y in k if int(y[1:]) == 1] # col 1
        bg_keys = bg_keys + [y for y in k if (int(y[1:]) == 12 and y not in bg_keys)]
        bg_keys = bg_keys + [y for y in k if (y[0] == 'A' and y not in bg_keys)]
        bg_keys = bg_keys + [y for y in k if (y[0] == 'H' and y not in bg_keys)]

        return bg_keys

    def get_data_keys(self,df):
        """Gets the dataframe keys for the data assuming a 1-well moat

        Args:
            df (pandas dataframe): datafram containing raw OD data

        Returns:
            list: list of keys
        """
        bg_keys = self.get_background_keys(df)

        data_keys = [k for k in df.keys() if k not in bg_keys]
        data_keys = data_keys[2:]

        return data_keys

    def estimate_background(self,df):
        """Estimates the OD background assuming a 1-well moat

        Args:
            df (pandas dataframe): datafram containing raw OD data

        Returns:
            float: background OD
        """
        bg_keys = self.get_background_keys(df)
        s = 0

        for key in bg_keys:
            s += np.average(df[key])

        bg = s/len(bg_keys)

        return bg


    def subtract_background(self,df):
        """Subtracts OD background from raw OD data

        Args:
            df (pandas dataframe): datafram containing raw OD data

        Returns:
            pandas dataframe: background-subtracted data
        """
        bg = self.estimate_background(df)
        datakeys = df.keys()[2:]

        if (df.values < bg).any():
            bg = np.min(df.values)
            

        for key in datakeys:
            df[key] = df[key] - bg

        return df

    def get_growth_rate_data(self,data_path):
        """Loads and background subtracts growth rate data

        Args:
            data_path (str): path to raw csv data

        Returns:
            pandas dataframe: background-subtracted raw data
        """
        df = dir_manager.load_growth_rate_data(data_path)
        df = self.subtract_background(df)

        return df

    def get_growth_rates_from_df(self,df,carrying_cap=None):
        
        """Estimates the growth rates from timeseries growth data in a dataframe

        Arguments:
            df: pandas DataFrame
                data frame containing growth rate data

        Returns:
            growth_rates: dict
                dictionary of growth rates for each experimental condition
        """

        growth_rates = {}

        # assumes outer row is background data
        data_keys = self.get_data_keys(df)
        time = df['Time [s]']

        for k in data_keys:
            growth_rates[k] = self.est_growth_rate(df[k],t=time)
            # cur_genotype = k[0]
            # if cur_genotype == prev_genotype:
            #     growth_rates[k] = est_growth_rate(df[k],t=time,carrying_cap=cc)
            # else:
            #     # get new carrying capacity
            #     genotype_keys = [k for k in data_keys if k[0]==prev_genotype]


        return growth_rates

    def gen_seascape_library(self,pop=None,debug=False):
        """Fits raw estimated growth rate values to a Hill dose-response curve

        Args:
            pop (population class object, optional): population class object. Defaults to None.
            debug (bool, optional): generates plots useful for debugging if true. Defaults to False.

        Raises:
            ValueError: raises error if there is no growth rate library in the population class

        Returns:
            dict: seascape library
        """
        if pop is None:
            pop = self

        if not 'growth_rate_library' in pop.__dict__:
            raise ValueError('No growth rate library in population.')
        else:
            gl = pop.growth_rate_library
        
        sl = {}
        dc = gl['drug_conc']
        
        for g in range(pop.n_genotype):
            popt = self.fit_hill_curve(dc,gl[str(g)],debug=debug)
            ic50 = popt[0]
            g_drugless = popt[1]
            hill_coeff = popt[2]
            d_t = {'ic50':ic50,
                'g_drugless':g_drugless,
                'hill_coeff':hill_coeff}
            sl[str(g)] = d_t

        return sl
        
    def sl_to_fitness(self,g,conc,pop=None):
        """Seascape library to fitness value (in units per second)

        Args:
            pop (population class object): population class object
            g (int): genotype
            conc (float): drug concentration

        Returns:
            float: fitness
        """

        if pop is None:
            pop = self

        ic50 = pop.seascape_library[str(g)]['ic50']
        
        g_drugless = pop.seascape_library[str(g)]['g_drugless']
        hc = pop.seascape_library[str(g)]['hill_coeff']

        f = self.logistic_pharm_curve(conc,ic50,g_drugless,hc)
        return f

    def gen_growth_rate_library(self,pop=None):

        if pop is None:
            # if type(self) is Population:
            #     pop = self
            # else:
            #     raise TypeError('Population object required')
            pop = self


        growth_rate_lib = {}
        growth_rate_lib['drug_conc'] = pop.seascape_drug_conc # add drug concentration data

        library = ['B','C','D','E','F'] # library of column names
        genotype = 0

        for df in pop.growth_rate_data: # for each plate
            
            growth_rates = self.get_growth_rates_from_df(df,carrying_cap=pop.max_od) # estimate the growth rate for each set of timeseries data

            for l in library:
                gr_vect = []
                i = 2
                for c in pop.seascape_drug_conc:
                    key = l + str(i) # each key is a seperate column in the plate reader data,looping over i holds genotype steady and increments drug concentration
                    gr_vect.append(growth_rates[key])
                    i+=1
                # print(len(gr_vect))
                growth_rate_lib[str(genotype)] = gr_vect
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
        
        growth_rate_lib[str(genotype)] = gr_vect

        return growth_rate_lib

    def est_growth_rate(self,growth_curve,t=None,debug=False,save_debug=False,num=None,
                        method='scipy',carrying_cap=4):

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
        growth_curve = growth_curve/carrying_cap

        if method == 'custom':
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
            
            r = np.max(r)

        # curve fit using scipy
        elif method == 'scipy':

            p0 = [10**-6,0.05,1]

            popt, pcov = scipy.optimize.curve_fit(self.logistic_growth_curve,
                                                t,growth_curve,p0=p0,
                                                bounds=(0,1))

            r = popt[0]
            if r < 0:
                r = 0
            if popt[2] < popt[1]: # if the carrying capacity is less than the initial population size
                r = 0
            if popt[2] < 0.05:
                r = 0
            
            if debug:
                fig,ax = plt.subplots()

                ax.plot(t,growth_curve)

                est = self.logistic_growth_curve(t,popt[0],popt[1],popt[2])
                
                ax.plot(t,est)
                # print(popt[0])
                p0 = round(popt[1]*10**5)/10**5
                k = round(popt[2]*10**5)/10**5
                r = round(r*10**5)/10**5
                title = str(r*(60**2)) + ' ' + str(p0) + ' ' + str(k)
                ax.set_title(title)
        
        if save_debug:
            if num is None:
                savename = 'debug' + os.sep + 'debug_fig.pdf'
            else:
                savename = 'debug' + os.sep + 'debug_fig' + str(num) + '.pdf'
            results_manager.save_fig(fig,savename)
        

        return r

    def fit_hill_curve(self,xdata,ydata,debug=False):

        # interpolate data
        xd_t = xdata
        yd_t = ydata
        f = scipy.interpolate.interp1d(xdata,ydata)

        if min(xdata) == 0:
            xmin = np.log10(xdata[1])
        else:
            xmin = np.log10(min(xdata))
        xmax = np.log10(max(xdata))

        xdata = np.logspace(xmin,xmax)
        if not xdata[0] == 0:
            xdata = np.insert(xdata,0,0)

        ydata = f(xdata)

        
        p0 = [0,ydata[0],-0.08]

        if ydata[0] == 0:
            g_drugless_bound = [0,1]
        else:
            # want the estimated drugless growth rate to be very close to the value given in ydata
            g_drugless_bound = [ydata[0]-0.0001*ydata[0],ydata[0]+0.0001*ydata[0]]

        bounds = ([-5,g_drugless_bound[0],-1],[4,g_drugless_bound[1],-0.001])
        popt, pcov = scipy.optimize.curve_fit(self.logistic_pharm_curve_vectorized,
                                            xdata,ydata,p0=p0,bounds=bounds)

        if debug:
            est = [self.logistic_pharm_curve(x,popt[0],popt[1],popt[2]) for x in xdata]
            fig,ax = plt.subplots()
            ax.plot(xd_t,yd_t)
            ax.plot(xdata,ydata)
            ax.plot(xdata,est)
            ax.set_xscale('log')
            ax.set_title('IC50 = ' + str(popt[0]))

        return popt

    def get_max_od(self,pop=None):
        """Extract the max OD from a plate

        Args:
            pop (population class, optional): _description_. Defaults to None.

        Returns:
            _type_: _description_
        """
        pop = self
            
        max_od = 0
        for df in pop.growth_rate_data:
            data_keys = self.get_data_keys(df)
            for key in data_keys:
                if max_od < max(df[key]):
                    max_od = max(df[key])

        return max_od

    def logistic_growth_curve(self,t,r,p0,k):
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

    def logistic_pharm_curve(self,x,IC50,g_drugless,hill_coeff):
        """Logistic dose-response curve. use if input is a single drug concentration

        Args:
            x (float): drug concentration scalar
            IC50 (float)): IC50
            g_drugless (float): drugless growth rate
            hill_coeff (float): Hill coefficient

        Returns:
            numpy array: array of growth rates
        """
        if x == 0:
            g = g_drugless
        else:
            g = g_drugless/(1+np.exp((IC50-np.log10(x))/hill_coeff))

        return g

    def logistic_pharm_curve_vectorized(self,x,IC50,g_drugless,hill_coeff):
        """Defines the logistic dose-response curve. Use if the input is a vector of drug concentration curves

        Args:
            x (numpy array): drug concentration vector
            IC50 (float)): IC50
            g_drugless (float): drugless growth rate
            hill_coeff (float): Hill coefficient

        Returns:
            numpy array: array of growth rates
        """
        g = []

        for x_t in x:
            if x_t == 0:
                g.append(g_drugless)
            else:
                g.append(g_drugless/(1+np.exp((IC50-np.log10(x_t))/hill_coeff)))

        return g