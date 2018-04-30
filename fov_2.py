import numpy as np
from plotting import *
from fov import K, N_fov, R_fov

def Mlim(Mhost, Mv1, Mv2):
    return Mhost*(10**((Mv2-Mv1)/-2.5))

print K(0.76, 1.5)

colors = ['royalblue', 'darkorange', 'mediumseagreen', 'crimson', 'darkorchid', 'turquoise']
Ds = [794., 880., 968.]

for D in Ds:

    # test for HSC at one M33 distance
    fov = 1.5

    if D == 794.:
        hNfields = [10.,23., 50., 100., 130.]
        mNfields = [10., 52., 100., 200., 293.]
    if D == 880.:
        hNfields = [10.,23., 50., 100., 106.]
        mNfields = [10., 52., 100., 200., 238.]
    if D == 968.:
        hNfields = [10.,23., 50., 88.]
        mNfields = [10., 52., 100., 197.]

    #M_v translates to Mstar translates to Nlum from LCDM
    rfov = R_fov(fov, D)
    rvir = 156.
    print D
    print 'max fields to observe ALL predicted sats at this D:', (rvir*0.76/rfov)**2.
    nlums = [10.5107, 10.1041, 9.6503, 9.15155, 8.60745, 8.03555, 7.42715, 6.7938, 6.1579, 5.51975, 4.89235, 4.28215, 3.71745, 3.186, 2.71315, 2.29685, 1.9339,  1.6118,  1.3349]#, 1.1006, 0.90715, 0.74345, 0.60935, 0.49755, 0.40475] #[10.5, 7.4, 3.7, 1.3]

    nlo = [ 6.99334539, 6.65277516, 6.28038182, 5.88619147, 5.45470996, 4.985979, 4.50332579, 4.00594141, 3.50884118, 3.01815135, 2.55360312, 2.11005008, 1.69224528, 1.32447856, 1.00788004, 0.72873106, 0.49636247, 0.30851167, 0.15671181]#, 0.03627285 -0.05957562 -0.13237366 -0.17996786 -0.21153889 -0.23183063]

    nhi = [ 14.06075461, 13.56412484, 13.03211818, 12.44580853, 11.79029004, 11.080321 , 10.33237421,  9.55585859,  8.75205882,  7.95434865,  7.16719688,  6.40544992,  5.67945472,  5.01192144,  4.39851996,  3.83386894,  3.33763753,  2.89468833,  2.50068819]#,  2.15612715,  1.85457562,  1.59877366,  1.38206786,  1.18913889,  1.02633063] 

    Mvs2 = [np.log10(1e3/3.2e9)*-2.5 + -18.8, np.log10(1e4/3.2e9)*-2.5 + -18.8, np.log10(1e5/3.2e9)*-2.5 + -18.8, np.log10(1e6/3.2e9)*-2.5 + -18.8]
    mstars =[1.00000000e+03,   1.46779927e+03,   2.15443469e+03,   3.16227766e+03, 4.64158883e+03,   6.81292069e+03,   1.00000000e+04,   1.46779927e+04, 2.15443469e+04,   3.16227766e+04,   4.64158883e+04,   6.81292069e+04,1.00000000e+05,   1.46779927e+05,   2.15443469e+05,   3.16227766e+05, 4.64158883e+05,   6.81292069e+05,   1.00000000e+06]#,   1.46779927e+06, 2.15443469e+06,   3.16227766e+06,   4.64158883e+06,   6.81292069e+06, 1.00000000e+07]  #[1e3, 1e4, 1e5, 1e6]
    Mvs = [np.log10(m/3.2e9)*-2.5 + -18.8 for m in mstars]
    #print Mvs, Mvs2
    
    plt.figure()
    ax1 = plt.subplot(111)
    ax2 = ax1.twiny()
    ax1.plot(Mvs, nlums, label=r'$\Lambda$CDM simulation', lw=5, color='gray')
    ax1.fill_between(Mvs, nlo, nhi, color='gray', alpha=0.25)

    obs = []
    for Mv,nlum in zip( Mvs,nlums):
        obs.append(nlum*K(rfov/rvir, 1.5))
    print 1,obs    
    ax1.plot(Mvs, obs,  c='k', label='1 pointing')
    ax1.legend()
    plt.savefig('M33_nfov_HSC_test.pdf')

    labels = ['', '(PAndAS)', '', '', r'($\rm R_{vir}$)']
    for i,c,l in zip(hNfields, colors, labels):
        obs = []
        print 'R_proj:', np.sqrt((i*rfov**2.)/rvir**2.)*rvir
        for Mv,nlum in zip(Mvs,nlums):
            obs.append(nlum*K(np.sqrt((i*rfov**2.)/rvir**2.), 1.5))
        print i, [round(o, 2) for o in obs]
        ax1.plot(Mvs, obs,  color=c, label='%i fields %s'%(i,l))

    ax1.legend()


    # add MegaCam
    rfov = R_fov(1., D)
    print rfov
    print D
    print 'max fields to observe ALL predicted sats at this D:', (rvir*0.76/rfov)**2.
    obs = []
    for Mv,nlum in zip(Mvs,nlums):
        obs.append(nlum*K(rfov/rvir, 1.5))
    print 1, obs
    #plt.plot(Mvs, obs, label='1 field', c='gray', ls='--')


    for i,c,l in zip(mNfields, colors, labels):
        obs = []
        print 'R_proj:', np.sqrt((i*rfov**2.)/rvir**2.)*rvir
        for Mv,nlum in zip(Mvs,nlums):
            obs.append(nlum*K(np.sqrt((i*rfov**2.)/rvir**2.), 1.5))
        print i, [round(o, 2) for o in obs]
        ax1.plot(Mvs, obs, label='%i fields %s'%(i,l), color=c, marker='o', ms=3,ls='--')
    ax1.legend(frameon=False)




    ax2.set_xticks(Mvs2)
    ax2.set_xlabel(r'limiting $\rm M_{*} \, (M_{\odot})$', fontsize=16)
    ax2.set_xlim(-10.5, -2.)
    ax2.set_xticklabels([r"$10^{3}$",r"$10^{4}$",r"$10^{5}$",r"$10^{6}$"])

    ax1.set_xlabel(r'limiting magnitude ($\rm M_V$)', fontsize=16)
    ax1.set_ylabel(r'total number of observed satellites', fontsize=16)
    ax1.axvline(x=-6.5, ymin=0., ymax=0.75,color='orange', zorder=-100, ls=':')
    ax1.set_xlim(-10.5, -2.)

    plt.figtext(0.05, 0.95, r'$\rm D_{M33} = %i \, kpc$'%D, color='black')
    plt.figtext(0.45, 0.8, 'solid: HSC', fontsize=14)
    plt.figtext(0.45, 0.75, 'dashed: MegaCam', fontsize=14)

    plt.savefig('M33_summary_sat_predictions_%i.pdf'%D)


