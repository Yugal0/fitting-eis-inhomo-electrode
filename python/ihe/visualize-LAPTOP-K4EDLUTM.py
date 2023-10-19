import matplotlib.pyplot as plt
import numpy as np

def plot_fit_profile(fitted_circuit,ax):
    ax=ax
    fig, ax=plt.subplots()
    
    if(fitted_circuit.name=="TLMtwo"):
        ax.plot([0,delta_0,delta_0,1],[t1_0,t1_0,t2_0,t2_0],linewidth=linewidth_actual,color=color_actual,label="Actual")
        ax.plot([0,delta_f1,delta_f1,1],[t1_f1,t1_f1,t2_f1,t2_f1],linewidth=linewidth,color=color_fitted1,label="Fitted 1")

    if(fitted_circuit.name=="TLMthree"):
        ax.plot([0,delta1_0,delta1_0,delta2_0+delta1_0,delta2_0+delta1_0,1],[t1_0,t1_0,t2_0,t2_0,t3_0,t3_0],linewidth=linewidth_actual,color=color_actual,label="Actual")
        ax.plot([0,delta1_f1,delta1_f1,delta2_f1+delta1_f1,delta2_f1+delta1_f1,1],[t1_f1,t1_f1,t2_f1,t2_f1,t3_f1,t3_f1],linewidth=linewidth,color=color_fitted1,label="Fitted 1")
    
    if(fitted_circuit.name=="TLMlin" or fitted_circuit.name=="TLMlinzert"):
        ax.plot([0,1],[t1_f1,t2_f1],linewidth=linewidth,color=color_fitted1,label="Fit 1")
        ax.plot([0,1],[t1_f2,t2_f2],linewidth=linewidth,color=color_fitted2,label="Fit 2")
    
    
    if(fitted_circuit.name=="TLMilin" or fitted_circuit.name=="TLMilinzert"):
        def give_ilin_R(r1,r2,delta):
            return 1/(1/r1+(1/r2-1/r1)*delta)
        dlt=np.linspace(0,1,50)
        ax.plot(dlt,give_ilin_R(t1_f1,t2_f1,dlt),linewidth=linewidth,color=color_fitted1,label="Fitted 1")
        ax.plot(dlt,give_ilin_R(t1_f2,t2_f2,dlt),linewidth=linewidth,color=color_fitted2,label="Fitted 2")

    
def plot_fit_error(z_in,z_fit,freq,ax):
    

def plot_fit_all(z_in,z_fit,freq,ax):