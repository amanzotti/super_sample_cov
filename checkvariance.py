import numpy as np
import matplotlib.pyplot as plt
import scipy.special
import scipy.interpolate

array=np.genfromtxt('/Users/alessandromanzotti/Downloads/camb/test_scal2Cls.dat', dtype = float)



ell=array[:,0]
Cphi=array[:,4]/(7.4311e12*array[:,0]**4)


thetaprecision=(0.0006/1.041)**2 ### theta is 1.041+0.0006


areatot=np.arange(1,10000)
sumwindow1=np.zeros(10000)
sumwindow2=np.zeros(10000)
sumapprox=np.zeros(10000)

for x in xrange(1,9800):
        theta=np.sqrt(areatot[x]/np.pi)*(np.pi/180)
        tck =scipy.interpolate.splrep(ell[:],ell[:]**2*(ell[:]+1)**2/4*(Cphi[:])*(ell[:]/(2*np.pi))*(2*scipy.special.j1(ell[:]*theta)/(ell[:]*theta))**2)
        sumwindow1[x]=scipy.interpolate.splint(1,2000,tck)
	sumwindow2[x]=np.sum((2*ell[:]+1)/(4*np.pi)*ell[:]**2*(ell[:]+1)**2/4*(Cphi[:])*4*(scipy.special.j1(ell[:]*theta)/(ell[:]*theta))**2)
        lmax=np.ceil(2*np.pi/theta)
        sumapprox[x]=np.sum((2*ell[1:lmax]+1)/(4*np.pi)*ell[1:lmax]**2*(ell[1:lmax]+1)**2/4*(Cphi[1:lmax]))
	#print np.sqrt(areatot[x]/np.pi)*(np.pi/180)
	#print lmax



fg = plt.figure(figsize=(10,8))

plt.rc('text', usetex=True)
plt.rc('font', family='serif')



plt.loglog(ell[:],ell[:]**2*(ell[:]+1)**2/4*(Cphi[:]),ell[:],scipy.interpolate.splev(ell[:], tck),'--')
plt.xlabel(r'$\ell$',fontsize=19)
plt.ylabel(r'$C^{\kappa}=\ell^{2}(\ell+1)^{2}C^{\phi}/4  $(K$^{2}$)',fontsize=17)
plt.ylim([5*10**-8,5*10**-7])
plt.xlim([2,200])
plt.savefig('Ckappa.pdf')
plt.clf()

cumsum2=np.cumsum((2*ell[:]+1)/(4*np.pi)*ell[:]**2*(ell[:]+1)**2/4*(Cphi[:]))
#plt.loglog(4*np.pi*180**2/(ell[:]**2),cumsum2
plt.loglog(4*np.pi*180**2/(ell[:]**2),thetaprecision*ell[:]/ell[:],label='Planck accuracy')
plt.loglog(np.arange(1,10001),sumwindow2[:],label='Integral')
plt.loglog(np.arange(1,10001),sumapprox[:],label='No window-l$_{max}$')
plt.loglog(np.arange(1,10001),sumwindow1[:],label='Sum approx')
plt.xlabel(r'Survey Area deg$^{2}$',fontsize=18)
plt.ylabel(r'$\sigma^{\kappa}_{lmax}$',fontsize=18)

legend = plt.legend(loc='best')
plt.ylim([6*10^-8,10**-3])
#plt2 = plt.twiny()
#plt2.set_xticklabels(tick_function(plt.xticks()))
#plt2.set_xlabel(r"Modified x-axis: $1/(1+X)$")
plt.savefig('variance.pdf')

plt.clf()

theta=np.sqrt(1000/np.pi)*(np.pi/180)
plt.loglog((ell[:100]),4*(scipy.special.j1(ell[:100]*theta)/(ell[:100]*theta))**2,label='Area 1000 $deg^{2}$')

theta=np.sqrt(600/np.pi)*(np.pi/180)
plt.loglog((ell[:100]),4*(scipy.special.j1(ell[:100]*theta)/(ell[:100]*theta))**2,label='Area 600 $deg^{2}$')

theta=np.sqrt(800/np.pi)*(np.pi/180)
plt.loglog((ell[:100]),4*(scipy.special.j1(ell[:100]*theta)/(ell[:100]*theta))**2,label='Area 800 $deg^{2}$')

plt.xlabel(r'$\ell$',fontsize=18)
plt.ylabel(r'$|\tilde W(l)|^2$',fontsize=19)
plt.xlim([1,100])
plt.ylim([6*10^-3,2])

legend = plt.legend(loc='best')


plt.savefig('window.pdf')




#plt.show()





