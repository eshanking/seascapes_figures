import numpy as np
import matplotlib.pyplot as plt

def hill(conc,ic50,g_drugless,hc):
    g = g_drugless/(1+np.exp((np.log10(ic50)-np.log10(conc))/hc))
    return g

def logistic_growth(rate,t,K,p0):
    L = (K-p0)/p0
    p = K/(1+(L)*np.exp(-rate*t))
    return p

fig,ax = plt.subplots(figsize=(4,4))

time = np.linspace(0,15,num=100)

for rate in [0.5,0.75,1,1.25]:
    pop = []
    for t in time:
        pop.append(logistic_growth(rate,t,1,0.01))
    ax.plot(time,pop,label=str(rate))

vline_y = np.linspace(0,1,num=10)
vline_x = np.ones(10)*12.5 
ax.plot(vline_x,vline_y,'--',color='black',linewidth=2)

vline_y = np.linspace(0,1,num=10)
vline_x = np.ones(10)*6 
ax.plot(vline_x,vline_y,'--',color='black',linewidth=2)

ax.legend(frameon=False)

ax.set_ylabel('Normalized OD')
ax.set_xlabel('Time (hr)')

fig.savefig('single_od_cartoon.pdf',bbox_inches='tight')