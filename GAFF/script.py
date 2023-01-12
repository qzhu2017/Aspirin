# This python scripts computes the overlap between two distributions
import numpy as np
np.set_printoptions(precision=4)

data1=np.genfromtxt("ACSALA/histo1")
data2=np.genfromtxt("ACSALA/histo2")
data3=np.genfromtxt("ACSALA13/histo1")
data4=np.genfromtxt("ACSALA13/histo2")
print('{:.4f}'.format(np.trapz(np.minimum(data1[:,1],data3[:,1]),x=data1[:,0])))
print('{:.4f}'.format(np.trapz(np.minimum(data2[:,1],data4[:,1]),x=data2[:,0])))

