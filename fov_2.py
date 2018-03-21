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
    nlums = [10.5, 7.4, 3.7, 1.3]
    Mvs = [np.log10(1e3/3.2e9)*-2.5 + -18.8, np.log10(1e4/3.2e9)*-2.5 + -18.8, np.log10(1e5/3.2e9)*-2.5 + -18.8, np.log10(1e6/3.2e9)*-2.5 + -18.8]
    mstars = [1e3, 1e4, 1e5, 1e6]


    plt.figure()
    ax1 = plt.subplot(111)
    ax2 = ax1.twiny()
    plt.plot(Mvs, nlums, label='LCDM simulation', lw=4, color='gray')

    obs = []
    for Mv,nlum in zip( Mvs,nlums):
        obs.append(nlum*K(rfov/rvir, 1.5))
    print 1,obs
    
    plt.plot(Mvs, obs,  c='k', label='1 pointing')

    plt.legend()
    plt.savefig('M33_nfov_HSC_test.pdf')

    labels = ['', '(PAndAS)', '', '', r'($\rm R_{vir}$)']
    for i,c,l in zip(hNfields, colors, labels):
        obs = []
        print 'R_proj:', np.sqrt((i*rfov**2.)/rvir**2.)*rvir
        for Mv,nlum in zip(Mvs,nlums):
            obs.append(nlum*K(np.sqrt((i*rfov**2.)/rvir**2.), 1.5))
        print i, [round(o, 2) for o in obs]
        plt.plot(Mvs, obs,  color=c, label='%i fields %s'%(i,l))

    plt.legend()
    plt.savefig('M33_nfov_HSC_test.pdf')

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
        plt.plot(Mvs, obs, label='%i fields %s'%(i,l), color=c, ls='--')
    plt.legend(frameon=False)
    ax1.set_xlabel(r'limiting magnitude ($\rm M_V$)', fontsize=16)
    ax1.set_xlim(-10, -2.5)
    ax2.set_xticks(Mvs)
    ax2.set_xticklabels([r"$10^{3}$",r"$10^{4}$",r"$10^{5}$",r"$10^{6}$"])
    ax2.set_xlabel(r'limiting $\rm M_{*} \, (M_{\odot})$', fontsize=16)
    ax1.set_ylabel(r'total number of observed satellites', fontsize=16)
    plt.figtext(0.05, 0.95, r'$\rm D_{M33} = %i \, kpc$'%D, color='black')
    plt.figtext(0.45, 0.8, 'solid: HSC', fontsize=14)
    plt.figtext(0.45, 0.75, 'dashed: MegaCam', fontsize=14)
    plt.axvline(x=-6.5, ymin=0., ymax=0.75,color='orange', zorder=-100, ls=':')
    plt.savefig('M33_summary_sat_predictions_%i.pdf'%D)


