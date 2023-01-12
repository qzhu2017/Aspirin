import numpy as np
import matplotlib.pyplot as plt

sigma=np.arange(0.04,0.08,0.01)
overlap=np.zeros(sigma.shape[0])
plt.rcParams["figure.figsize"] = (3.5*2,7*2)
fig, axs= plt.subplots(sigma.shape[0],sharex=True) #,sharey=True)

for i in range(sigma.shape[0]):
    string='{:.2f}'.format(sigma[i])
    d1=np.genfromtxt("ACSALA/histo1-"+string)
    d2=np.genfromtxt("ACSALA13/histo1-"+string)
    axs[i].plot(d1[:,0], d1[:,1])
    axs[i].plot(d2[:,0], d2[:,1])
    axs[i].set_xlim([-0.4,2])
    axs[i].text(-0.35, np.amax(d2[:,1])/2,r"$\sigma$=" + string,fontsize=15)
    #print(1.8-i*1.9/sigma.shape[0],20.0-i*20.0/sigma.shape[0])
#axs[5].set_ylabel("Probability density")
axs[sigma.shape[0]-1].set_xlabel(r"$k(\chi,\chi')$")

#plt.tight_layout()
plt.savefig('dist1.png')

fig, axs= plt.subplots(sigma.shape[0],sharex=True) #,sharey=True)
for i in range(sigma.shape[0]):
    string='{:.2f}'.format(sigma[i])
    d3=np.genfromtxt("ACSALA/histo2-"+string)
    d4=np.genfromtxt("ACSALA13/histo2-"+string)
    axs[i].plot(d3[:,0], d3[:,1])
    axs[i].plot(d4[:,0], d4[:,1])
    axs[i].set_xlim([-0.4,2])
    axs[i].text(-0.35, np.amax(d4[:,1])/2,r"$\sigma$=" + string,fontsize=15)
axs[sigma.shape[0]-1].set_xlabel(r"$k(\chi,\chi')$")
plt.savefig('dist2.png')
