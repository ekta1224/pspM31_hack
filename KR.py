import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

Z = 1.5
rvir = 101.
r = 33.
R = r/rvir


def K(R):
    k1 = -0.2615
    k2 = 6.888
    k3 = -7.035
    k4 = 0.9667
    k5 = 0.5298
    k6 = 0.2055

    if R < 0.2:
        #print 'R < 0.2'
        return (k1*R + k2*R**2. + k3*R**3.)

    if R >= 0.2:
        #print 'R > 0.2'
        return (k4*np.arctan((R/k5) - k6))
        

def K2(Z, R, A=True, B=False):
    k1 = -0.2615
    k2 = 6.888
    k3 = -7.035
    k4 = 0.9667
    k5 = 0.5298
    k6 = 0.2055
    if A:
        #print 'using method A'
        if R < 0.2:
            return 0.5 * (R**2. + Z**2.)**(-0.5) * (k1 + 2.*k2*R + 3.*k3*R**2.) 
        #return (R)/(R**2 + Z**2)**(2.) * (k1 + 2.*k2*R + 3.*k3*R**2.) 
        #return 1./(R*(R**2 + Z**2)) * (k1 + 2.*k2*R + 3.*k3*R**2.) 
        if R >= 0.2:
            return 0.5 * (R**2. + Z**2.)**(-0.5) * (k4/(k5*((k6 - (R/k4))**2.+ 1)))
        #return (R)/(R**2 + Z**2)**(2.) * (k4/(k5*((k6 - (R/k4))**2.+ 1)))
        #return 1./(R*(R**2 + Z**2)) * (k4/(k5*((k6 - (R/k4))**2.+ 1)))

    if B:
        #print 'using method B'
        if R < 0.2:
            return 0.5 * (1./R) * (R**2. + Z**2.)**(-1.) * (k1 + 2.*k2*R + 3.*k3*R**2.) * (1./rvir**2.)
        #return (R)/(R**2 + Z**2)**(2.) * (k1 + 2.*k2*R + 3.*k3*R**2.) 
        #return 1./(R*(R**2 + Z**2)) * (k1 + 2.*k2*R + 3.*k3*R**2.) 
        if R >= 0.2:
            return 0.5 * (1./R) * (R**2. + Z**2.)**(-1.) * (k4/(k5*((k6 - (R/k4))**2.+ 1))) * (1./rvir**2.)


print 'R and K(R)', R, K(R)

Ks = []
K2s = []
for R in np.arange(0.01, 1.5, .05):
    myK = integrate.dblquad(K2, 0., R, lambda Z: 0., lambda Z: 1.)[0]
    myK2 = integrate.dblquad(K2, 0., R, lambda Z: 0., lambda Z:1.5)[0]
    #print R, myK, myK2
    Ks.append(myK)
    K2s.append(myK2)

plt.figure()
ax = plt.subplot(111)
plt.plot(np.arange(0.01, 1.5, .05), Ks, label='Z=1')
plt.plot(np.arange(0.01, 1.5, .05), K2s, label='Z=1.5')
plt.plot(0.33, 0.54,'o ', color='orange')
plt.plot(0.5, 0.76,'bo')
ax.grid(color='k', linestyle=':', linewidth=1)
plt.legend()
plt.savefig('R_vs_Klos_A_halfZ.pdf')
