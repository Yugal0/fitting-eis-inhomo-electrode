import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp

def give_z_for_arbitary_tau(tau,freq_list,ax,name):
    
    def fun(x,y,lmd):
        return np.vstack((y[2]*tau(x),y[3]*tau(x),-lmd*y[1],lmd*y[0]))

    def bc(ya,yb):
        return np.array([ya[0]-1,yb[2],yb[3],ya[1]])

    Q=6.61e-4
    Rref=194.56

    x=np.linspace(0,1,1000)
    y=np.zeros((4,x.size))

    Z_list=np.zeros((np.size(freq_list)),dtype=complex)
    ii=0
    for f in freq_list:
        w=2*np.pi*f
        lmd2=w*Q*Rref
        sol1 = solve_bvp(lambda x,y: fun(x,y,lmd=lmd2), bc, x, y)
        y3=sol1.sol(x)[2]
        y4=sol1.sol(x)[3]
        Z=2/(-(1/(Rref))*(y3[0]+y4[0]*1j)) 
        Z_list[ii]=Z
        ii=ii+1

    Z_list_real=np.real(Z_list)
    Z_list_imag=np.imag(Z_list)

    ax.plot(Z_list_real/(Rref/1.5),-Z_list_imag/(Rref/1.5),label=name)
    ax.set_aspect("equal")
    ax.set_xlim([0,6])
    ax.set_ylim([0,10])
    
    return Z_list

def tau_plotter(tau,ax,name):
    x=np.linspace(0,1,100)
    ax.plot(x,tau(x),label=name)