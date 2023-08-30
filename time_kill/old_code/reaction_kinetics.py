import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import scipy

def kinetic_eqn(y,t,ka,kd,N):
    # print(y)
    r,R = y
    ka = 10**ka
    kd = 10**kd
    N = 10**N
    dydt = [-ka*r*N,ka*r*N - kd*R]
    return dydt

ka = -8
kd = -5

y0 = [1,0]

t = np.linspace(0,80000,num=1000)

# fig,ax = plt.subplots()

N_vect = [3,4,5,6,7,8]

ydata = []

for N in N_vect:
    sol = odeint(kinetic_eqn,y0,t,args=(ka,kd,N))
    ydata.append(sol[:,1])
    # ax.plot(t,sol[:,1],label='%.1E' % N)

# ax.legend(frameon=False,title='Population size',loc='upper right')
# ax.set_ylabel('Normalized intensity')
# ax.set_xlabel('Time (s)')

def kinetic_sol(t,y0,ka,kd,N):
    y = odeint(kinetic_eqn,y0,t,args=(ka,kd,N))
    return y

# Define a global error function that sums up the squared errors for each data set 
def global_loss(params,ydata,t):

    ka,kd,Ns = params[0],params[1],params[2:]
    indx = 0
    loss = 0

    for yd in ydata:
        y_t = kinetic_sol(t,[1,0],ka,kd,Ns[indx])[:,1]
        loss += (np.array(y_t - yd)**2).sum()
        indx+=1
    
    return loss

N_guess = np.ones(len(N_vect))*6
# N_guess = [4,7]
x0 = np.concatenate(([-8],[-5],N_guess))

res = scipy.optimize.minimize(global_loss,x0,args=(ydata,t),method='Nelder-Mead',
                              options={'disp':True,'adaptive':True,'xatol':0.01})

fig,ax_list = plt.subplots(nrows=len(N_vect))
ka = res.x[0]
kd = res.x[1]

N_est = res.x[2:]

indx = 0
for N in N_est:
    ax_list[indx].plot(t,ydata[indx],color='black')
    y_t = kinetic_sol(t,y0,ka,kd,N)[:,1]
    ax_list[indx].plot(t,y_t,color='r')
    ax_list[indx].annotate(str(N),xy=(50000,0))
    ax_list[indx].set_ylim(-0.1,1.1)
    indx +=1
