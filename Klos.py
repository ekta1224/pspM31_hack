import numpy as np
from scipy.integrate import dblquad
from scipy import interpolate


def get_numerical_normed_los(R,Z):
    # This is K(R)  using the raw values which were used in creaing the fitting
    # function of Eq. 4. Beyond R=1, the values are an extrapolation, and become increasingly less accurate.
    # It was computed on the discrete values r_over_Rvir.
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
    tck =  interpolate.splrep(r_over_Rvir, K_of_R, k=1)

    def dkdr(r):
        return interpolate.splev(r,tck,der=1)

    t0 = np.arctan2(R,Z)
    term1 = dblquad(lambda r, theta: dkdr(r)*np.sin(theta), 0, t0, lambda x: 0, lambda x: Z/np.cos(x))
    term2 = dblquad(lambda r, theta: dkdr(r)*np.sin(theta), t0, np.pi/2, lambda x: 0, lambda x: R/np.sin(x))
    print 'finished an iteration of integrating', R
    print term1[0] + term2[0], 'is K(R, Z) for R=', R, 'Z=', Z
    return term1[0] + term2[0]


get_numerical_normed_los(.4, 1.5)

# Answer to check against
#Z=1.5
#R = [0.0, 0.05,.1,.2, .3,.4,.5,.6,.7,.8,.9, 1, 1.1, 1.2, 1.3, 1.4,1.5]
#K(R,Z) = [0.0, 0.02928724, 0.10360976, 0.30315803, 0.495490586338, 0.64919431, 0.770664268383, 0.86710225, 0.946403823717, 1.00873837, 1.06091190553, 1.103422, 1.14177189992, 1.17847986, 1.21349950337, 1.24713824, 1.27914254]

