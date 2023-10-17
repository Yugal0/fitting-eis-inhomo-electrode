import models
from impedance.models.circuits import CustomCircuit
import numpy as np

c=CustomCircuit(initial_guess=[6000,0.6,0.2,0.2,0.200407010357445,0.7,1],circuit="TLMthree")
freq=np.logspace(-3,3,100)
cc=c.predict(freq)
print(cc)
