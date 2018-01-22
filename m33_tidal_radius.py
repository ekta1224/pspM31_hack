import numpy as np
import matplotlib.pyplot as plt
from plotting import *
from cosmo_tools import r_vir
from nfw import NFW

rvir =  r_vir(0.3, 0.7, 1e11)

Rs = np.arange(0., rvir*3., 20.)

mratios = [1., 1/2.,1/5.,  1/10., 1/50., 1/100., 1/200.]

plt.figure()
ax = plt.subplot(111)
for m in mratios:
    rts = []
    for r in Rs:
        rts.append(r*m**(1./3.)) #m=sat/host mass
    plt.plot(Rs, rts, label='mass ratio %s'%m)
#ax.tick_params(axis='x',which='minor',bottom='on')
#ax.tick_params(axis='y',which='minor',bottom='on')
plt.minorticks_on()
plt.xlabel('host-satellite separation [kpc]')
plt.ylabel('M33 tidal radius')
plt.title('for Mvir=1e11 Msun, Rvir~%s kpc'%int(rvir))
plt.xlim(0,200)
plt.ylim(0,200)
plt.legend(loc='best')
plt.savefig('m33_rt.png')


#tidal radius of M31
rvir = 330.
Rs = np.arange(0., rvir*2., 20.)
m = 0.1

plt.figure()
plt.plot(Rs, [r*m**(1./3) for r in Rs])
plt.xlabel('host-satellite separation [kpc]')
plt.ylabel('M31 tidal radius')
plt.title('for Mvir=1.5e12 Msun, Rvir~%s kpc'%int(rvir))
plt.xlim(0,400)
plt.ylim(0,200)

plt.savefig('m31_rt.png')



#update after the M33 hack day satellite project - 9/18/2017
m33mass = 1e11
m31mass = 1.5e12
m31c = 9.56

def rt(peri):
    m31 = NFW(m31mass, peri, m31c)
    m31enc = m31.mass()
    return peri*(2.*m31enc/m33mass)**(-1./3.), m31enc

rs = np.arange(10., 300., 5.)

masses = []
rts = []
for r in rs:
    thisrt = rt(r)
    rts.append(thisrt[0])
    masses.append(thisrt[1])

print rt(50.), rt(100.), rt(150.)

plt.figure()
plt.scatter(rts, rs, c=masses, edgecolor='black', norm=matplotlib.colors.LogNorm())
plt.ylabel(r'$\rm r_{peri, M33}$ [kpc]', fontsize=15)
plt.xlabel(r'$\rm r_{tidal, M33}$ [kpc]', fontsize=15)
plt.axhline(y=55, color='lime', label='McConnachie+09 pericenter')
#ax1 = plt.gca()
#ax2 = ax1.twinx()
#ax2.set_yscale("log")
#ax2.plot(rts,masses)
cbar = plt.colorbar()
cbar2 = cbar.set_label(r'$\rm M_{M31}(r_{peri, M33})$ [M$_{\odot}$]', size=15)
#cbar.tick_params(labelsize=16) 
plt.legend()
plt.savefig('m33_rt_updated.pdf')
