from impedance.models.circuits.elements import element
import numpy as np
import scipy
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