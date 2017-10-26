import numpy as np
import matplotlib.pyplot as plt
import scipy.special
import scipy.interpolate


## load lenisng potential from CAMB
array=np.genfromtxt('/home/manzotti/Software/camb/test_scalCls.dat', dtype = float)


## remember its normalization and it is also multiplied by a l^4 


ell=array[:,0]
Cphi=array[:,4]/(7.4311e12*array[:,0]**4)


thetaprecision=(0.0006/1.041)**2 ### theta is 1.041+0.0006


areatot=np.arange(1,10000)
sumwindow1=np.zeros(10000)
sumwindow2=np.zeros(10000)
sumapprox=np.zeros(10000)


## compute sigmak of mode outside the box. It is done via a sum or an integral. The result does not change that much but check yourself.

for x in xrange(1,9800):
        theta=np.sqrt(areatot[x]/np.pi)*(np.pi/180)
        fsky=(areatot[x]*(np.pi/180)**2)/(4*np.pi)
        tck =scipy.interpolate.splrep(ell[:],ell[:]**2*(ell[:]+1)**2/4*(Cphi[:])*(ell[:]/(2*np.pi))*(2*scipy.special.j1(ell[:]*theta)/(ell[:]*theta))**2)
        sumwindow1[x]=fsky*scipy.interpolate.splint(1,4000,tck)
	sumwindow2[x]=fsky*np.sum((2*ell[:]+1)/(4*np.pi)*ell[:]**2*(ell[:]+1)**2/4*(Cphi[:])*4*(scipy.special.j1(ell[:]*theta)/(ell[:]*theta))**2)
        lmax=np.ceil(2*np.pi/theta)
        sumapprox[x]=fsky*np.sum((2*ell[1:lmax]+1)/(4*np.pi)*ell[1:lmax]**2*(ell[1:lmax]+1)**2/4*(Cphi[1:lmax]))
	#print np.sqrt(areatot[x]/np.pi)*(np.pi/180)
	#print lmax


print areatot[40]

varcom=np.genfromtxt('/media/USB20FD/Wayne-B_lensing/Code/variancenonoisenofsky_fisher.txt', dtype = float)
#varcomno=np.genfromtxt('/media/USB20FD/Wayne-B_lensing/Code/variance3000nonoise_fisher.txt', dtype = float)

fg = plt.figure(figsize=(10,8))
ax1 = fg.add_subplot(111)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


# or sp = plt.subplot2grid((3,2),(0,0))
# l1=plt.loglog(array[2:80,0],(array[2:80,3])**2/(array[2:80,2]),array[2:80,0],(array[2:80,4]),'--')

#p1=plt.plot(array[:,0],abs(W1[:])/(0.02*5000*(60**2)),array[:,0],array[:,1]/(0.02*5000*(60**2)),array[:,0],abs(W2[:])/(0.02*5000*(60**2)),array[:,0],abs(W3[:])/(0.02*5000*(60**2)),array[:,0],abs(W4[:])/(0.02*5000*(60**2)))

ax1.plot(varcom[:,0],(varcom[:,1]),varcom[:,0],(varcom[:,2]),varcom[:,0],(varcom[:,3]))
plt.xlabel(r'$\ell_{\rm{max}}$',fontsize=18)
plt.ylabel(r'$ \sigma^{2}_{s}\times f_{sky}$ ',fontsize=18)
plt.tick_params(axis='both', which='major', labelsize=18)
plt.tick_params(axis='both', which='minor', labelsize=18)
#plt.xlim((4,3000))
#ax1.set_xscale('log')
ax1.set_yscale('log')
#plt.text(1, 1, r'$f_{\rm{sky}}=0.024$', fontsize=15)
#plt.title( r'$f_{\rm{sky}}=0.024$')


ax1.plot(np.arange(1,10001),np.arange(1,10001)/np.arange(1,10001)*sumwindow1[100],label=r'$\sigma_{\kappa}^{2} \times f_{sky}~~$100 deg2')
ax1.plot(np.arange(1,10001),np.arange(1,10001)/np.arange(1,10001)*sumwindow1[700],label=r'$\sigma_{\kappa}^{2} \times f_{sky}~~$700 deg2')
ax1.plot(np.arange(1,10001),np.arange(1,10001)/np.arange(1,10001)*sumwindow1[1000],label=r'$\sigma_{\kappa}^{2} \times f_{sky}~~$1000 deg2')
ax1.plot(np.arange(1,10001),np.arange(1,10001)/np.arange(1,10001)*sumwindow1[2000],label=r'$\sigma_{\kappa}^{2} \times f_{sky}~~$2000 deg2')

legend = ax1.legend()
plt.ylim(10**(-9),10**(-6))
plt.xlim(800,4000)

plt.savefig('variancecompare.pdf')
plt.clf()
fg = plt.figure(figsize=(10,8))
ax1 = fg.add_subplot(111)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
print sumwindow1[600]
ax1.plot(np.arange(1,10001),sumwindow1[:],100,5.95*10**-6,'ro',600,6.72*10**-7,'ro',1000,4.73*10**-7,'ro',2500,9.5*10**-8,'ro')
plt.xlabel(r'$Sky Area (sq .deg)$',fontsize=18)
plt.ylabel(r'$ \sigma_{\kappa}^{2}$ ',fontsize=18)
ax1.set_xscale('log')
ax1.set_yscale('log')
plt.savefig('prova.pdf')
