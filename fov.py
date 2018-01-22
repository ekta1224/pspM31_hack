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

def R_fov(fov, distance):
    return fov*np.pi/180.*distance/2.

print N_fov(7.4,  R_fov(1.5,800)/115.,1.5) #satellites per pointing
