import numpy as np
import matplotlib.pyplot as plt
import scipy.special

array=np.genfromtxt('/media/USB20FD/Wayne-B_lensing/Code/derivatives2.txt', dtype = float)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
adj = plt.subplots_adjust(hspace=0.7,wspace=0.7)

# or sp = plt.subplot2grid((3,2),(0,0))
# l1=plt.loglog(array[2:80,0],(array[2:80,3])**2/(array[2:80,2]),array[2:80,0],(array[2:80,4]),'--')

#p1=plt.plot(array[:,0],abs(W1[:])/(0.02*5000*(60**2)),array[:,0],array[:,1]/(0.02*5000*(60**2)),array[:,0],abs(W2[:])/(0.02*5000*(60**2)),array[:,0],abs(W3[:])/(0.02*5000*(60**2)),array[:,0],abs(W4[:])/(0.02*5000*(60**2)))
print max(array[:,3])
p1=plt.plot(array[:,0],(array[:,1]),'r',array[:,0],array[:,3],'b',array[:,0],array[:,5],'k')
plt.xlabel(r'$\ell$',fontsize=14)
plt.ylabel(r'$\frac{\partial C_{\ell}}{\partial s}=C^{\rm fid}\frac{\partial \ln (l^2 C_l^{\rm fid}) }{\partial \ln l}$ ',fontsize=16)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.tick_params(axis='both', which='minor', labelsize=20)
#plt.ylim((0,50000))
plt.xlim((0,3000))

plt.savefig('derivatives.pdf')
