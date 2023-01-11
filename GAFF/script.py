# This python scripts computes the overlap between two distributions
import numpy as np
data1=np.genfromtxt("ACSALA/histo")
data2=np.genfromtxt("ACSALA13/histo")
print(np.trapz(np.minimum(data1[:,1],data2[:,1]),x=data1[:,0]))

