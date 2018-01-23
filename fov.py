# compute the multiplicative constant for FOV calculations
# Dooley 2017a, Equation 6
# number of satellites along the line of sight relative to the number expected within Rvir 
import numpy as np

def K(R,Z):
    k1 = -0.0440
    k2 = 0.3913
    k3 = 0.9965
    k4 = -0.3438
    
    return (k2*Z/2.)*np.log(1. +(R**2./Z**2.)) + (k3*Z**2.)*(np.sqrt(1+(R**2./Z**2.)) - 1.) + (3.*k4*R**2.*Z/2.) + k2*R*(np.pi/2. - np.arctan(R/Z)) + k3*R**2.*np.log(Z/R*(1. + np.sqrt(1. + (R**2./Z**2.)))) + k1

def N_fov(Nlum, R, Z):
    return K(R,Z)*Nlum

def R_fov(fov, distance, sq=False):
    if sq:
        return fov *(np.pi/180.)**2.*distance/2.
    return fov*np.pi/180.*distance/2.

m33r = 960. 
rvir = 120.
Z = 1.5
nlum = 7.4
fov = 2.2

print 'K(R):', K(R_fov(fov, m33r )/rvir, Z)
#print N_fov(nlum, R_fov(fov, m33r)/rvir, Z) #satellites per pointing

fovs = [1., 1.5, 2.2]
nlums = [10.5, 7.4, 3.7, 1.3]
rs = [794., 880., 960.]
for nlum in nlums:
    print N_fov(nlum, R_fov(fov, m33r)/rvir, Z) #satellites per pointing

#below the R_fov of hypersuprimecam (R~0.06), K pretty much goes to negative numbers
rvir = 120.
Z = 1.5
nlum = 10.5
fov = 2.2

for m33r in rs:
    print m33r
    for fov in fovs:
        print fov
        rfov = R_fov(fov, m33r)
        k = K(R_fov(fov, m33r)/rvir, Z)
        nfov = N_fov(nlum, rfov/rvir,Z)
        points = rvir**2./R_fov(fov, m33r)**2.
        print ' K, R_fov,  Npoint %s & %s & %s & %s & %s & %s & %s'%(round(k, 3), round(rfov,2), int(round(points,0)), round(N_fov(10.5, rfov/rvir, Z),3), round(N_fov(7.4, rfov/rvir, Z),3), round(N_fov(3.7, rfov/rvir, Z),3), round(N_fov(1.3, rfov/rvir, Z),3))
        minfov = 120*0.06*2.*180/np.pi/m33r
        #print 'minimum FOV (deg) as fn of m33r', minfov
    
print N_fov(10.5, R_fov(0.85, 960.)/rvir, Z) #satellites per pointing


