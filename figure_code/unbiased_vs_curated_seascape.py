import numpy as np
import matplotlib.pyplot as plt

x = np.arange(-10,10,0.1)

y = (x-10)**2*(x-1)**2*(x+4)*(x+10)
# y = (x+5)**2*(x+1)**2*(x-7)

fig,ax = plt.subplots()
ax.plot(x,y)
ax.set_xlim(-5,10)
ax.set_ylim(-1000,70000)