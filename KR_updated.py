import numpy as np
from scipy import integrate
from scipy import interpolate
import matplotlib.pyplot as plt

# Z = 1.5
# rvir = 101.
# r = 33.
# R = r/rvir

r_over_Rvir = np.arange(0,2.001,.05)
K_of_R = np.array([ 0.        ,  0.004329  ,  0.03463203,  0.09307359,  0.16693723,
        0.24837662,  0.33495671,  0.41477273,  0.4862013 ,  0.55357143,
        0.6155303 ,  0.66856061,  0.71888528,  0.76704545,  0.81385281,
        0.85227273,  0.88555195,  0.92099567,  0.95183983,  0.9767316 ,
        1.        ,  1.02      ,  1.04      ,  1.06      ,  1.08      ,
        1.1       ,  1.12      ,  1.14      ,  1.16      ,  1.18      ,
        1.2       ,  1.22      ,  1.24      ,  1.26      ,  1.28      ,
        1.3       ,  1.32      ,  1.34      ,  1.36      ,  1.38      ,
        1.4       ])
print len(r_over_Rvir), len(K_of_R)

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

thisKR = [K(r) for r in r_over_Rvir]

plt.figure()
plt.plot(r_over_Rvir, K_of_R, label='GD')
plt.plot(r_over_Rvir, thisKR, label='EP')
plt.legend()
plt.savefig('KR_comparison.pdf')


def mydKdr(R):
    k1 = -0.2615
    k2 = 6.888
    k3 = -7.035
    k4 = 0.9667
    k5 = 0.5298
    k6 = 0.2055

    if R < 0.2:
        return k1 + 2.*k2*R + 3*k3*R**2.

    if R >= 0.2:
        return k4 / (k5 * ((k6 - (R/k4))**2. + 1))

def myKlos(R,Z):
    theta0 = np.arctan2(R,Z)
    term1 = integrate.dblquad(lambda r, theta: mydKdr(r)*np.sin(theta), 0, theta0, lambda x: 0., lambda x: Z/np.cos(x))
    term2 = integrate.dblquad(lambda r, theta: mydKdr(r)*np.sin(theta), theta0, np.pi/2., lambda x: 0., lambda x: R/np.sin(x))
    return term1[0] + term2[0]


def Klos(R,Z):
    Rs = np.arange(0., 2., .05)
    y =  [K(r) for r in Rs]
    spl = interpolate.splrep(Rs, y, k=1)

    def dKdr(R):
        fn = interpolate.splev(R, spl, der=1)
        return fn
    
    theta0 = np.arctan2(R,Z)
    term1 = integrate.dblquad(lambda r, theta: dKdr(r)*np.sin(theta), 0, theta0, lambda x: 0., lambda x: Z/np.cos(x))
    term2 = integrate.dblquad(lambda r, theta: dKdr(r)*np.sin(theta), theta0, np.pi/2., lambda x: 0., lambda x: R/np.sin(x))
    return term1[0] + term2[0]

                              
Ks = []
K2s = []
myKs = []
myK2s = []
Rs = np.arange(0., 2., .05)
for R in Rs:
    myK = Klos(R, Z=1.)
    myK2 = Klos(R, Z=1.5)
    print R, myK, myK2
    Ks.append(myK)
    K2s.append(myK2)
    myKs.append(myKlos(R, Z=1.))
    myK2s.append(myKlos(R, Z=1.5))
                
plt.figure()
ax = plt.subplot(111)
plt.plot(Rs, Ks, label='Z=1')
plt.plot(Rs, K2s, label='Z=1.5')
plt.plot(Rs, myKs, label='Z=1 (analytic)')
plt.plot(Rs, myK2s, label='Z=1.5 (analytic')
plt.plot(0.33, 0.54,'o ', color='orange')
plt.plot(0.5, 0.76,'bo')
ax.grid(color='k', linestyle=':', linewidth=1)
plt.legend()
plt.xlim(0,1.5)
plt.savefig('Klos_updated.pdf')
