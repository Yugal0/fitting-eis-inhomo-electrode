import numpy as np
def find_dip_angle(fitted_circuit):
    p=fitted_circuit.parameters_
    
    dip_angle=alpha*45-min(np.angle(z))