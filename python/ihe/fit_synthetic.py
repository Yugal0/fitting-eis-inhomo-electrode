import numpy as np

def fit_synthetic(circuit_simu,params_sim,circuit_fit,initial_gs_fit):



def add_gaussian_noise(z,err):
    #err is the relative SD from mean value
    z_with_noise=np.random.normal(np.real(z),err*np.abs(np.real(z)))+1j*np.random.normal(np.imag(z),err*np.abs(np.imag(z)))
    return z_with_noise
