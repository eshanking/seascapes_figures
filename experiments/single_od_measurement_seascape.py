from fears.utils import AutoRate
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as sciopt

def hill(conc,ic50,g_drugless,hc):
    g = []
    for c in conc:
        if c <= 0:
            g.append(g_drugless)
        else:
            g.append(g_drugless/(1+np.exp((np.log10(ic50)-np.log10(c))/hc)))
    return g

def hill_scalar(conc,ic50,g_drugless,hc):
    if conc <= 0:
        g = g_drugless
    else:
        g = g_drugless/(1+np.exp((np.log10(ic50)-np.log10(conc))/hc))
    return g

folder_path = '/Users/eshanking/repos/fears/fears/data/od_plates_no_lid'

layout_path = '/Users/eshanking/repos/fears/fears/data/plate_layout.csv'
reference_path = '/Users/eshanking/repos/fears/fears/data/plates/20210929_plate1.csv'

t_obs = 12*3600
drug_conc = [0, 0.01, 1, 10, 100, 1000]

e = AutoRate.Experiment(folder_path=folder_path,mode='single_measurement',
                        exp_layout_path=layout_path,
                        ref_data_path=reference_path,t_obs=t_obs,drug_conc=drug_conc)
e.execute()

fig = plt.figure(figsize=(8,8))
plot_counter = 1

for g in e.growth_rate_lib.keys():
    grd = e.growth_rate_lib[g]
    gr_vect = []
    gr_err = []
    for c in grd.keys():
        gr_vect.append(grd[c]['avg'])
        gr_err.append(grd[c]['std'])

    gr_vect = [g*3600 for g in gr_vect]
    gr_err = [g*3600 for g in gr_err]

    # remove drug-free condition
    gr_vect = gr_vect[1:]
    gr_err = gr_err[1:]

    gr_vect = gr_vect/max(gr_vect)

    p0 = [0,0.99,-1]
    # bounds = ([-5,5],[0.9,1.1],[-2,0])
    bounds = ([-10**5,0.9,-3],[10**5,1.1,0])
    popt, pcov = sciopt.curve_fit(hill,
                                drug_conc[1:],gr_vect,
                                p0=p0,
                                bounds=bounds)

    ic50 = popt[0]
    g_drugless = popt[1]
    hc = popt[2]

    conc_t = np.logspace(-3.1,3.1,100)
    fit = [hill_scalar(c,ic50,g_drugless,hc) for c in conc_t]

    plt.subplot(4,4,plot_counter)
    # plt.scatter(drug_conc,gr_vect)
    plt.errorbar(drug_conc[1:],gr_vect,yerr=gr_err,fmt='o')
    plt.plot(conc_t,fit)
    plt.scatter(ic50,hill_scalar(ic50,ic50,g_drugless,hc),color='r')
    plt.xscale('log')
    # plt.ylim(0,0.5)
    plt.xlim(10**-3,10**5)
    plt.title(str(round(hc,2)),fontsize=15)
    plot_counter+=1


axes = fig.axes
axes[0].set_ylabel('OD',fontsize=15)
axes[4].set_ylabel('OD',fontsize=15)
axes[8].set_ylabel('OD',fontsize=15)
axes[12].set_ylabel('OD',fontsize=15)

for a in axes:
    a.tick_params(axis='both', which='major', labelsize=15)

fig.tight_layout()
fig.savefig('single_od_seascape.pdf')