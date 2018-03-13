import numpy as np
from plotting import *
from fov import K, N_fov, R_fov

def Mlim(Mhost, Mv1, Mv2):
    return Mhost*(10**((Mv2-Mv1)/-2.5))


colors = ['royalblue', 'darkorange', 'mediumseagreen', 'crimson', 'darkorchid', 'turquoise']

# test for HSC at one M33 distance
fov = 1.5
D = 794.

#M_v translates to Mstar translates to Nlum from LCDM
rfov = R_fov(fov, D)
print rfov
rvir = 156.
nlums = [10.5, 7.4, 3.7, 1.3]
Mvs = [np.log10(1e3/3.2e9)*-2.5 + -18.8, np.log10(1e4/3.2e9)*-2.5 + -18.8, np.log10(1e5/3.2e9)*-2.5 + -18.8, np.log10(1e6/3.2e9)*-2.5 + -18.8]
mstars = [1e3, 1e4, 1e5, 1e6]


plt.figure()
ax1 = plt.subplot(111)
ax2 = ax1.twiny()
plt.plot(Mvs, nlums, label='LCDM simulation', lw=3, color='k')

obs = []
for Mv,nlum in zip(Mvs,nlums):
    obs.append(nlum*K(rfov/rvir, 1.5))
plt.plot(Mvs, obs,  c='gray') #label='1 pointing',

plt.legend()
plt.savefig('M33_nfov_HSC_test.pdf')

for i,c in zip([10., 50., 100., 150.],colors):
    obs = []
    for Mv,nlum in zip(Mvs,nlums):
        obs.append(nlum*K(np.sqrt((i*rfov**2.)/rvir**2.), 1.5))
    plt.plot(Mvs, obs,  color=c)#label='%i pointings'%(i),

plt.legend()
plt.savefig('M33_nfov_HSC_test.pdf')

# add MegaCam
rfov = R_fov(1., D)
print rfov
obs = []
for Mv,nlum in zip(Mvs,nlums):
    obs.append(nlum*K(rfov/rvir, 1.5))
plt.plot(Mvs, obs, label='1 pointing', c='gray', ls='--')

for i,c in zip([10., 50., 100., 150., 200., 250.], colors):
    obs = []
    for Mv,nlum in zip(Mvs,nlums):
        obs.append(nlum*K(np.sqrt((i*rfov**2.)/rvir**2.), 1.5))
    plt.plot(Mvs, obs, label='%i pointings'%(i), ls='--', color=c)

plt.legend()
ax1.set_xlabel(r'limiting magnitude ($\rm M_V$)', fontsize=16)
ax1.set_xlim(-10, -2.5)
ax2.set_xticks(Mvs)
ax2.set_xticklabels([r"$10^{3}$",r"$10^{4}$",r"$10^{5}$",r"$10^{6}$"])
ax2.set_xlabel(r'limiting $\rm M_{*} \, (M_{\odot})$', fontsize=16)
ax1.set_ylabel(r'total number of observed satellites', fontsize=16)
plt.figtext(0.05, 0.95, r'$\rm D_{M33} = %i \, kpc$'%D, color='red')
plt.figtext(0.45, 0.8, 'solid: HSC', fontsize=14)
plt.figtext(0.45, 0.75, 'dashed: MegaCam', fontsize=14)
plt.savefig('M33_summary_sat_predictions_794.pdf')


