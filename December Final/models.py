from impedance.models.circuits.elements import element
from impedance.models.circuits import CustomCircuit
import numpy as np
import scipy
import matplotlib.pyplot as plt

@element(num_params=2, units=["Ohm", "F"],overwrite=True)
def TLMuni(p, f):
   
    omega = 2 * np.pi * np.array(f)
    R1, Q = p[0], p[1]
    Z= 2*np.sqrt(R1/(omega*1j*Q))*(1/(np.tanh(np.sqrt(1j*omega*R1*Q))))
    return Z
@element(num_params=5, units=["Ohm", "Ohm F", "Ohm F","",""],overwrite=True)
def TLMtwo(p, f):
   
    omega = 2 * np.pi * np.array(f)
    R1, R1Q, R2Q, delta1, alpha = p[0], p[1], p[2], p[3], p[4]
    delta2=1-delta1
    k12=R2Q/R1Q
    wc1=1/(R1Q)
    wc2=1/(R2Q)
    S1=np.sqrt((omega*1j)**alpha/wc1)
    S2=np.sqrt((omega*1j)**alpha/wc2)
    C1=1/np.tanh(delta1*S1)
    C2=1/np.tanh(delta2*S2)

    Z=2*(R1)*(C1*C2*S1*k12+S2)/(S1*(C2*S1*k12+C1*S2))

    return Z

@element(num_params=7, units=["Ohm", "Ohm F", "Ohm F", "","", "",""],overwrite=True)
def TLMthree(p, f):
    omega = 2 * np.pi * np.array(f)
    R1, R1Q, R2Q, R3Q, delta1, delta2, alpha = p[0], p[1], p[2],p[3], p[4], p[5], p[6]
    delta3=1-delta1-delta2
    delta2=delta2
    Q=R1Q/R1
    R2=R2Q/Q
    R3=R3Q/Q
    k12=R2/R1
    k13=R3/R1
    wc1=1/(R1*Q)
    wc2=1/(R2*Q)
    wc3=1/(R3*Q)
    S1=np.sqrt((omega*1j)**alpha/wc1)
    S2=np.sqrt((omega*1j)**alpha/wc2)
    S3=np.sqrt((omega*1j)**alpha/wc3)
    C1=1/np.tanh(delta1*S1)
    C2=1/np.tanh(delta2*S2)
    C3=1/np.tanh(delta3*S3)

    Z=2*(R1)*((S1*S3*C1*k12**2+S2*C2*(S1*C1*C3*k13+S3)*k12+S2**2*C3*k13)/(S1*(S3*S1*k12**2+S2*C2*(S3*C1+S1*C3*k13)*k12+S2**2*C1*C3*k13)))

    return Z

@element(num_params=4, units=["Ohm", "Ohm F", "Ohm F", ""],overwrite=True)
def TLMlin(p,f):
    omega = 2 * np.pi * np.array(f)
    R1, R1Q, R2Q, alpha = p[0], p[1], p[2], p[3]
    Q=R1Q/R1
    R2=R2Q/Q
    t12=R1/R2
    S=np.sqrt((1j*omega)**alpha*R1*Q)
    kk1=-2/3*(S*(t12**(-1/2)/(t12-1)))
    kk2=-2/3*(t12*S/(t12-1))
    Z = (2*R1/S)*(scipy.special.iv(1/3, kk1)*scipy.special.iv(2/3, kk2)-scipy.special.iv(-1/3, kk1)*scipy.special.iv(-2/3, kk2))/((scipy.special.iv(-1/3, kk1))*scipy.special.iv(1/3, kk2)-scipy.special.iv(1/3, kk1)*scipy.special.iv(-1/3, kk2))
    return Z

@element(num_params=4, units=["Ohm", "Ohm F", "Ohm F", ""],overwrite=True)
def TLMlinzert(p,f):
    omega = 2 * np.pi * np.array(f)
    R1, R1Q, R2Q, alpha = p[0], p[1], p[2], p[3]
    Q=R1Q/R1
    R2=R2Q/Q
    eps=(R2-R1)/R1
    S=np.sqrt((R1*Q*(omega*1j)**alpha))

    #second expansion
    Z=2*4*R1*(np.exp(4*S)+2*np.exp(2*S)+1)*(np.exp(2*S)+1)*S/((((S**4-2/3*(S**3)+2*S**2-5/2*S+1/8)*eps**2+(4*S**3+S)*eps-4*S**2)*np.exp(2*S)+((-S**4-2/3*(S**3)-2*S**2-5/2*S-1/8)*eps**2+(4*S**3+S)*eps+4*S**2)*np.exp(4*S)+(4*S**2-S*eps+7/8*(eps**2))*np.exp(6*S)-4*S**2-S*eps-7*eps**2*(1/8)))
    return Z

@element(num_params=4, units=["Ohm", "Ohm F", "Ohm F", ""],overwrite=True)
def TLMilin(p,f):
    omega = 2 * np.pi * np.array(f)
    R1, R1Q, R2Q, alpha = p[0], p[1], p[2], p[3]
    Q=R1Q/R1
    R2=R2Q/Q
    t12=R1/R2
    wclin1=1/(R1*Q)
    wclin2=1/(R2*Q)
    S1=np.sqrt((1j*omega)**alpha/wclin1)
    S2=np.sqrt((1j*omega)**alpha/wclin2)
    
    Z = -(2*1j*R1/S1)*(((scipy.special.yv(1,-2*1j*t12*S2/(t12-1)))*(scipy.special.jv(0,-2*1j*S1/(t12-1))))-((scipy.special.jv(1,-2*1j*t12*S2/(t12-1)))*(scipy.special.yv(0,-2*1j*S1/(t12-1)))))/(((scipy.special.jv(1,-2*1j*t12*S2/(t12-1)))*(scipy.special.yv(1,-2*1j*S1/(t12-1))))-((scipy.special.yv(1,-2*1j*t12*S2/(t12-1)))*(scipy.special.jv(1,-2*1j*S1/(t12-1)))))
    
    return Z

@element(num_params=4, units=["Ohm", "Ohm F", "Ohm F", ""],overwrite=True)
def TLMilinzert(p,f):
    omega = 2 * np.pi * np.array(f)
    R1, R1Q, R2Q, alpha = p[0], p[1], p[2], p[3]
    Q=R1Q/R1
    R2=R2Q/Q
    wclin1=1/(R1*Q)
    eps=R1/R2-1
    S=np.sqrt((1j*omega)**alpha/wclin1)

    Z=2*R1*4*(np.exp(4*S)+2*np.exp(2*S)+1)*(np.exp(2*S)+1)*S/((((S**4+2*S**3+2*S**2+3/2*S+9/8)*eps**2+(-4*S**3-S)*eps-4*S**2)*np.exp(2*S)+((-S**4+2*S**3-2*S**2+3/2*S-9/8)*eps**2+(-4*S**3-S)*eps+4*S**2)*np.exp(4*S)+(-(1/8)*eps**2+S*eps+4*S**2)*np.exp(6*S)-4*S**2+S*eps+(1/8)*eps**2))
    
    return Z

def profile_plotter(circuit,ax):
    
    if(circuit._is_fit()):
        param=circuit.parameters_
    else:
        param=circuit.initial_guess
    profile=circuit.circuit

    ax1=ax
    if(profile=="TLMtwo"):
       R1=param[0]
       R1Q=param[1]
       R2Q=param[2]
       delta1=param[3]
       Q=R1Q/R1
       R2=R2Q/Q
       ax1.plot([0,delta1,delta1,1],[R1,R1,R2,R2],linewidth=2,color="red")
    
    if(profile=="TLMthree"):
       R1=param[0]
       R1Q=param[1]
       R2Q=param[2]
       R3Q=param[3]
       delta1=param[4]
       delta2=param[5]
       Q=R1Q/R1
       R2=R2Q/Q
       R3=R3Q/Q
       ax1.plot([0,delta1,delta1,delta2+delta1,delta2+delta1,1],[R1,R1,R2,R2,R3,R3],linewidth=2,color="blue")

    if(profile=="TLMlin" or profile=="TLMlinzert"):
       R1=param[0]
       R1Q=param[1]
       R2Q=param[2]
       R2=R1*R2Q/R1Q
       ax1.plot([0,1],[R1,R2],linewidth=2,color="green")

    if(profile=="TLMilin" or profile=="TLMilinzert"):
       R1=param[0]
       R1Q=param[1]
       R2Q=param[2]
       R2=R1*R2Q/R1Q
       ax1.plot([0,1],[R1,R2],linewidth=2,color="blue")
       def give_ilin_R(r1,r2,delta):
        return 1/(1/r1+(1/r2-1/r1)*delta)
       dlt=np.linspace(0,1,50)
       ax1.plot(dlt,give_ilin_R(R1,R2,dlt),linewidth=2,color="grey")


def error_plotter(freq,zdata1,zdata2,ax):
    err=np.abs(zdata1-zdata2)**2
    ax.semilogx(freq,err)
    ax.set_xlabel("freq")
    ax.set_ylabel("error")

def add_noise(z,err):
    #err is the relative SD from mean value
    z_with_noise=np.random.normal(np.real(z),err*np.abs(np.real(z)))+1j*np.random.normal(np.imag(z),err*np.abs(np.imag(z)))
    return z_with_noise

# def synthetic_checker(generater_model_name,gen_params=[],detector_model_name,detect_guess=[],freq,noise):
#     cgenerater=CustomCircuit(initial_guess=gen_params,circuit=generater_model_name)
#     zgen=cgenerater.predict(freq)
#     add_noise(zgen,noise)
#     cdetector=CustomCircuit(initial_guess=detect_guess,circuit=detector_model_name)

def give_weights(z,wt):
    if (wt=="mod"):
        return np.concatenate((np.abs(z),np.abs(z)))
    elif (wt=="prop"):
        return np.concatenate((np.real(z),np.imag(z)))
    elif (wt=="" or wt=="unit"):
        return np.ones(2*len(z))
    
def find_index_of_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx