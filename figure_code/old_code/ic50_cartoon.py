from seascapes_figures.classes.population_class import Population
from seascapes_figures.utils import results_manager
import matplotlib.pyplot as plt
import numpy as np

np.random.seed(2022)

p = Population()

conc = np.logspace(-4,4)

ydata = []
ic50 = 0
# convert to uM
ic50 = ic50-6

for c in conc:

    y_t = p.logistic_equation(c,1,ic50)
    ydata.append(y_t)

fig,ax = plt.subplots()

ax.plot(conc,ydata,color='black',linewidth=2)

# plot y intercept

xdata = np.logspace(-4,0)
ydata = 0.5*np.ones(len(xdata))

ax.plot(xdata,ydata,'--',color='gray')
# plot x intercept

ydata = np.arange(0,11)/20
xdata = np.ones(len(ydata))

ax.plot(xdata,ydata,'--',color='gray')

ax.set_xscale('log')

ax.set_ylim(0,1)
ax.set_xlim(10**-4,10**4)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

ax.set_ylabel('Growth rate (hr^-1)',fontsize=22)
ax.set_xlabel('Drug concentration (uM)',fontsize=22)

ax.annotate('IC50', xy=(2, 0.5),  xycoords='data',fontsize=18,
            xytext=(0.7, 0.6), textcoords='axes fraction',
            arrowprops=dict(facecolor='black', width=1))

ax.set_xticks([10**-3,10**-1,10**1,10**3])

plt.xticks(fontsize=22)
plt.yticks(fontsize=22)

results_manager.save_fig(fig,'ic50_cartoon.png')


xdata = np.random.normal(size=100)
ydata = []

for x in xdata:
    ydata.append(np.random.normal(loc=x,scale=0.7))

fig,ax = plt.subplots(figsize=(7,4))
ax.scatter(xdata,ydata,color='black')

ax.set_xlim([-3,3])
ax.set_ylim([-3,3])

ydata = xdata
ax.plot(xdata,ydata,color='red',label='PCA 1')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.legend(frameon=False,fontsize=20)

plt.xticks(fontsize=22)
plt.yticks(fontsize=22)

results_manager.save_fig(fig,'PCA_cartoon.png')