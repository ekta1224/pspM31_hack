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
        return (k1*R + k2*R**2. + k3*R**3.)/(4. * np.pi * R**2.)

    if R >= 0.2:
        #print 'R > 0.2'
        return (k4*np.arctan((R/k5) - k6))/(4. * np.pi * R**2.)
        
print 'R and K(R)', R, K(R)


print 'integration by parts'
#Ks = []
#for R in np.arange(0.1, 1.5, .1):

KK =  4.*np.pi*Z * ( ( K(R)/(4.*np.pi*R) ) + integrate.quad(K, 0., R)[0])
print R, KK, integrate.quad(K, 0., R)[0]
#Ks.append(KK)

#plt.figure()
#plt.plot(np.arange(0.1, 1.5, .1), Ks)
#plt.savefig('R_vs_Klos.pdf')
