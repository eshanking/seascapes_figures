import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# cmap = mpl.cm.get_cmap('copper')

# x = np.arange(-5,10,0.1)

# # y = (x-10)**2*(x-1)**2*(x+4)*(x+10)
# # y = (x+5)**2*(x+1)**2*(x-7)

# fig,ax = plt.subplots(figsize=(7,3.5))

# num = 1

# for i in [0,1,2,4]:

#     y = (x-10)**2*(x-i)**2*(x+5)

#     y = np.array(y)
#     y = y/np.max(y)

#     ax.plot(x,y,color=cmap(num/4),linewidth=3)
#     num += 1

# ax.set_xlim(-5,10)

# ax.tick_params(axis='y', labelsize=12)
# ax.set_xticks([])
# ax.set_xlabel('Genotype',fontsize=14)
# ax.set_ylabel('Fitness',fontsize=14)
# # ax.colorbar()

# cb = fig.colorbar(mpl.cm.ScalarMappable(norm=None, cmap=cmap), ax=ax,ticks=[])
# cb.set_label(label='env. gradient',fontsize=12)

# ax.set_ylim(0,10000)

#%% Surface plot

fig,ax = plt.subplots(figsize=(6,3))

# fig,ax = plt.subplots(figsize=(5,3))

g = np.arange(-5,10,0.1)
e = np.arange(0,5,0.1)

G, E = np.meshgrid(g,e)

F = (G-10)**2*(G-E)*(G+5)

F = F-np.min(F)
F = F/np.max(F)

im = ax.pcolormesh(G,E,F,cmap='magma')
cont = ax.contour(G,E,F,[0.4,0.8],colors='black')
ax.clabel(cont, fmt='%2.1f', colors='w', fontsize=14)
cb = fig.colorbar(im)
cb.set_label(label='Fitness',fontsize=14)

ax.set_yticks([])
ax.set_xticks([])
ax.set_xlabel('Genotype',fontsize=14)
ax.set_ylabel('Environment',fontsize=14)

fig.savefig('unbiased_vs_curated.pdf',bbox_inches='tight')
# %%
