import numpy as np
import matplotlib.pyplot as plt
import astropy.coordinates as coord
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from os.path import expanduser
home = expanduser("~")

#%%
def make_field_centers(sc, nrows):
    # 'sc' is an Astropy SkyCoord object with the central RA, Dec coordinate
    # 'nrows' is the number of rows of fields to be created in the Declination direction
    # 'nfields_width' is the number of fields in RA for the center row
    
    # NOTE: 'nrows' must be an odd number
    
# separation of field centers in degrees; HSC FOV ~1.5deg diameter; 10% field overlap
    overlap_frac_ra = 0.1
    offset_ra_deg = (1.5*(1.0-overlap_frac_ra))*u.deg 
    overlap_frac_dec = 0.2
    offset_dec_deg = (1.5*(1.0-overlap_frac_dec))*u.deg 
    
    ra_offset = offset_ra_deg/np.cos(sc.dec.radian)
    dec_offset = offset_dec_deg
    
    rownums = np.arange(nrows)-int(nrows/2)
    
    fields_ra = []
    fields_dec = []
    
# Row "0" is the one centered on the target. Odd-numbered rows should be shifted by a half-field offset
# The number of fields in each row is "nfields_width" for row 0, or generally: "nfields_width - N" for row N

    for row in rownums:
        if np.abs(row) >= 0:
            nfields_row = 15# np.abs(row+1)
            if np.mod(row, 2) == 0:
                rashifts = np.arange(nfields_row) - int(nfields_row/2)
            else:
                rashifts = np.arange(nfields_row) - int(nfields_row/2) + 0.5
            
            racen = sc.ra + (rashifts*ra_offset)
            deccen = racen*0.0 + (sc.dec - (row*dec_offset))
        
            fields_ra = np.append(fields_ra, racen.value)
            fields_dec = np.append(fields_dec, deccen.value)
    
    fields_ra = np.array(fields_ra)
    fields_dec = np.array(fields_dec)
  
    return fields_ra, fields_dec

#%%
def make_df_field_centers(sc, nrows):
    # 'sc' is an Astropy SkyCoord object with the central RA, Dec coordinate
    # 'nrows' is the number of rows of fields to be created in the Declination direction
    # 'nfields_width' is the number of fields in RA for the center row
    
    # NOTE: 'nrows' must be an odd number
    
# separation of field centers in degrees; DF FOV ~2.9x3.3; 10% field overlap
    overlap_frac_ra = 0.2
    offset_ra_deg = (2.9*(1.0-overlap_frac_ra))*u.deg 
    overlap_frac_dec = 0.1
    offset_dec_deg = (3.3*(1.0-overlap_frac_dec))*u.deg 
    
    ra_offset = offset_ra_deg/np.cos(sc.dec.radian)
    dec_offset = offset_dec_deg
    
    rownums = np.arange(nrows)-int(nrows/2)
    
    fields_ra = []
    fields_dec = []
    
# Row "0" is the one centered on the target. No shift applied to each row.
# The number of fields in each row is "nfields_width" for row 0, or generally: "nfields_width - N" for row N

    for row in rownums:
        if np.abs(row) >= 0:
            nfields_row = 11 # np.abs(row+1)
            rashifts = np.arange(nfields_row) - int(nfields_row/2)
            #if np.mod(row, 2) == 0:
            #    rashifts = np.arange(nfields_row) - int(nfields_row/2)
            #else:
            #    rashifts = np.arange(nfields_row) - int(nfields_row/2) + 0.5
            
            racen = sc.ra + (rashifts*ra_offset)
            deccen = racen*0.0 + (sc.dec - (row*dec_offset))
        
            fields_ra = np.append(fields_ra, racen.value)
            fields_dec = np.append(fields_dec, deccen.value)
    
    fields_ra = np.array(fields_ra)
    fields_dec = np.array(fields_dec)
  
    return fields_ra, fields_dec


#%%
# Suprime-Cam FOV: (~30x30') -- actually 34'x27'
def suprimecam_square(ra,dec):
    size=(30.0/2.0)/60.0
#    racorn = [(ra-size)/np.cos(dec*u.deg),(ra-size)/np.cos(dec*u.deg),(ra+size)/np.cos(dec*u.deg),(ra+size)/np.cos(dec*u.deg),(ra-size)/np.cos(dec*u.deg)]
    racorn = [(ra-size),(ra-size),(ra+size),(ra+size),(ra-size)]
    deccorn = [(dec-size),(dec+size),(dec+size),(dec-size),(dec-size)]
    return(racorn,deccorn)

#%%
# Dragonfly FOV: (~2.9x3.3 deg)
def dragonfly_fov(ra,dec):
    xsize = 2.9
    ysize = 3.3
#    racorn = [(ra-size)/np.cos(dec*u.deg),(ra-size)/np.cos(dec*u.deg),(ra+size)/np.cos(dec*u.deg),(ra+size)/np.cos(dec*u.deg),(ra-size)/np.cos(dec*u.deg)]
    racorn = [(ra-xsize/2),(ra-xsize/2),(ra+xsize/2),(ra+xsize/2),(ra-xsize/2)]
    deccorn = [(dec-ysize/2),(dec+ysize/2),(dec+ysize/2),(dec-ysize/2),(dec-ysize/2)]
    return(racorn,deccorn)

#%%
### How to use:

# 1. Uncomment the galaxy of interest below (both the skycoord and the name),
#    and comment out all the others.
# 2. Select the desired number of fields (in the call to 'make_field_centers'),
#    where the first input is the number of rows of fields, and the second is
#    the number of fields in the widest (center) row.
# 3. After the code runs, the target list will be output as 'targets.txt', and 
#    the fields plot will be saved as 'fields.png'. If you want these saved, 
#    rename them so they're not overwritten.

## M33:
sc0 = SkyCoord(ra='01h33m50.904s',dec='30d39m34.79s',distance=880.0*u.kpc,frame='icrs')
targetname = 'M33'


# To offset the center position so it's centered on a CCD:
#cen_offset = -0.3*u.deg # offset the center position by this amount to avoid chip gaps
#sc2dec = sc0.dec+cen_offset
#sc = SkyCoord(ra=sc0.ra,dec=sc2dec,distance=sc0.distance,frame='icrs')

# No offset applied:
sc = sc0

############################################
### Select the field centers:
ra_fields, dec_fields = make_field_centers(sc, 22)
#ra_fields, dec_fields = make_df_field_centers(sc, 9)

sc_fields_all = SkyCoord(ra = ra_fields*u.deg, dec = dec_fields*u.deg, frame='icrs')
sc_fields_sep = sc_fields_all.separation(sc0)

# Remove a few to fit within a 2-night observing request:
keep_fields = sc_fields_all.ra.value < 0 # create boolean array with all values True
#keep_fields[((sc_fields_sep.degree) < 3.5*170./50.)] = True
keep_fields[((sc_fields_sep.degree) < 3.5*170./50.) & ((sc_fields_sep.degree) > 2.0)] = True
#keep_fields[((sc_fields_sep.degree) < 3.5*100./50.) & ((sc_fields_sep.degree+0.75) > 3.5)] = True
pd_fields = (sc_fields_all.ra.value < 22.7) & (sc_fields_all.dec.value > 30.75)
keep_fields[pd_fields]=False

keep_fields20B = sc_fields_all.ra.value > 1000 # create boolean array with all values False
keep_fields20B[[186,187,201,202,215,216,230,231,245]] = True

#keep_fields19A = sc_fields_all.ra.value < 0 # create boolean array with all values True
#keep_fields19A[[186,187,201,202,216,217,231,232,233,246,247]] = True

sc_fields = sc_fields_all[keep_fields]
#sc_fields_not19A = sc_fields_all[keep_fields & ~keep_fields19A]
#sc_fields19A = sc_fields_all[keep_fields19A]
sc_fields_not20B = sc_fields_all[keep_fields & ~keep_fields20B]
sc_fields20B = sc_fields_all[keep_fields20B]


#%%

params = {
   'axes.labelsize': 24,
   'font.size': 16,
   'legend.fontsize': 10,
#   'xtick.labelsize': 16,
   'xtick.major.width': 3,
   'xtick.minor.width': 2,
   'xtick.major.size': 8,
   'xtick.minor.size': 5,
   'xtick.direction': 'in',
   'xtick.top': True,
   'lines.linewidth':3,
   'axes.linewidth':3,
   'axes.labelweight':3,
   'axes.titleweight':3,
   'ytick.major.width':3,
   'ytick.minor.width':2,
   'ytick.major.size': 8,
   'ytick.minor.size': 5,
   'ytick.direction': 'in',
   'ytick.right': True,
#   'ytick.labelsize': 20,
   'text.usetex': True,
   'text.latex.preamble': r'\boldmath',
   'figure.figsize': [14, 11]
   }

plt.rcParams.update(params)

from astropy.coordinates import CartesianRepresentation, SphericalRepresentation, Angle, SkyCoord

n = 10000
rvir = 161.
d = 794.
m33 = sc0

fieldcen = m33

fig = plt.figure(figsize=(8,8))
ax = plt.subplot(111)
plt.xlabel('RA [deg]', fontsize=18)
plt.ylabel('Dec [deg]', fontsize=18)
and22 = SkyCoord('01h27m40.0s','28d05m25s', unit='deg')
m33 = SkyCoord('01h33m50.0s','30d39m37s', unit='deg')
m31 = SkyCoord('0h42m44.3s', '41d16m9s', unit='deg')
plt.plot(and22.ra.deg, and22.dec.deg, '*', color='darkviolet' ,mec='white',ms=16, label='AND XXII', zorder=5)

#m.tissot(m33.ra.value,m33.dec.value,3.5,100,facecolor='none',edgecolor='k')
#m.tissot(m33.ra.value,m33.dec.value,3.5*100.0/50.0,100,facecolor='none',edgecolor='k',linewidth=8)

# unit circle at Dec=30 deg.
cpts = np.arange(0,360.0,1)
xpts = 1.0*np.cos(np.deg2rad(cpts))/np.cos(30.0*u.deg)
ypts = 1.0*np.sin(np.deg2rad(cpts))

xcirc50 = (xpts*3.5)+m33.ra.value
ycirc50 = (ypts*3.5)+m33.dec.value

plt.plot(xcirc50,ycirc50,color='cornflowerblue')

xcirc100 = (xpts*3.5*100.0/50.0)+m33.ra.value
ycirc100 = (ypts*3.5*100.0/50.0)+m33.dec.value

plt.plot(xcirc100,ycirc100,color='orange',ls='--')

xcirc161 = (xpts*3.5*161.0/50.0)+m33.ra.value
ycirc161 = (ypts*3.5*161.0/50.0)+m33.dec.value

plt.plot(xcirc161,ycirc161,color='green',ls=':')

dp = np.loadtxt('corners_PAndAS11.txt')
ra = dp[:,0]
dec = dp[:,1]
plt.plot(ra, dec, lw=2, color='k')
plt.xlim(*plt.xlim()[::-1])
plt.figtext(0.43, 0.7, 'PAndAS', fontsize=20, color='k')
circle= plt.Circle((m33.ra.deg, m33.dec.deg), 3.5, fill=False, color='cornflowerblue', lw=2, label='50 kpc')
circle2 = plt.Circle((m33.ra.deg, m33.dec.deg), 3.5*rvir/50., fill=False, color='green', lw=2, ls=':')
circle3 = plt.Circle((m33.ra.deg, m33.dec.deg), 3.5*100./50., fill=False, color='orange', lw=2, ls='--')

plt.xlim(39,-6)
plt.ylim(16,53)

#######
from matplotlib.patches import Circle, Wedge, Polygon, Ellipse, Rectangle
from matplotlib.collections import PatchCollection

ring = Wedge((m33.ra.deg, m33.dec.deg), 3.5*100./50., 0, 360, width=(3.5*100./50. - 3.5)),  # Full ring
p = PatchCollection(ring, alpha=0.1)
p.set_color('orange')

m31ell = Ellipse((m31.ra.deg, m31.dec.deg), (200./60.0), (50.0/60.0), angle=35, fill=True, linewidth=2, color='k')
ax.add_patch(m31ell)

m33ell = Ellipse((m33.ra.deg, m33.dec.deg), (73./60.0), (45.0/60.0), angle=23, fill=True, linewidth=2, color='k')
ax.add_patch(m33ell)

#m.draw_hsc_focal_planes(sc_fields_not20B.ra.value,sc_fields_not20B.dec.value,color='k',alpha=0.1)
#m.draw_hsc_focal_planes(sc_fields20B.ra.value,sc_fields20B.dec.value,color='b',alpha=0.4)
#m.draw_hsc_focal_planes(sc_fields_not19A.ra.value,sc_fields_not19A.dec.value,color='k',alpha=0.1)
#m.draw_hsc_focal_planes(sc_fields19A.ra.value,sc_fields19A.dec.value,color='b',alpha=0.4)

#rafield = 21.0
#decfield = 26.0
rasize = 2.9
decsize = 3.3
#rabox = [rafield-rasize/2, rafield-rasize/2, rafield+rasize/2, rafield+rasize/2, rafield-rasize/2]
#decbox = [decfield-decsize/2, decfield+decsize/2, decfield+decsize/2, decfield-decsize/2, decfield-decsize/2]

for i in range(len(sc_fields)):
    field = Rectangle([sc_fields[i].ra.deg-rasize/2, sc_fields[i].dec.deg-decsize/2], 2.9, 3.3, color='b', fill=False)
    ax.add_patch(field)

plt.legend(numpoints=1)
plt.figtext(0.77, 0.23, r'50 kpc', color='cornflowerblue')
plt.figtext(0.77, 0.2, r'100 kpc', color='orange')
plt.figtext(0.77, 0.17, r'161 kpc', color='green')

#plt.savefig('M33_HSC_fields.pdf')
plt.savefig('M33_HSC_fields_ALL_21B.pdf')

#%%
fig = plt.figure(figsize=(8,8))
ax = plt.subplot(111)
plt.xlabel('RA [deg]', fontsize=18)
plt.ylabel('Dec [deg]', fontsize=18)
plt.plot(and22.ra.deg, and22.dec.deg, '*', color='darkviolet' ,mec='white',ms=16, label='AND XXII', zorder=5)

plt.plot(ra, dec, lw=2, color='k')
plt.xlim(*plt.xlim()[::-1])
plt.figtext(0.77, 0.7, 'PAndAS', fontsize=20, color='k')
circle= plt.Circle((m33.ra.deg, m33.dec.deg), 3.5, fill=False, color='cornflowerblue', lw=2, label='50 kpc')
circle2 = plt.Circle((m33.ra.deg, m33.dec.deg), 3.5*rvir/50., fill=False, color='limegreen', lw=2, ls=':')
circle3 = plt.Circle((m33.ra.deg, m33.dec.deg), 3.5*100./50., fill=False, color='orange', lw=2, ls='--')

plt.plot(xcirc50,ycirc50,color='cornflowerblue')
plt.plot(xcirc100,ycirc100,color='orange',ls='--')

plt.xlim(29,17)
plt.ylim(23,34)

ring = Wedge((m33.ra.deg, m33.dec.deg), 3.5*100./50., 0, 360, width=(3.5*100./50. - 3.5)),  # Full ring
p = PatchCollection(ring, alpha=0.1)
p.set_color('orange')

m33ell = Ellipse((m33.ra.deg, m33.dec.deg), (73./60.0), (45.0/60.0), angle=23, fill=True, linewidth=2, color='k')
ax.add_patch(m33ell)

rafield = 21.0
decfield = 26.0
rasize = 2.9
decsize = 3.3
rabox = [rafield-rasize/2, rafield-rasize/2, rafield+rasize/2, rafield+rasize/2, rafield-rasize/2]
decbox = [decfield-decsize/2, decfield+decsize/2, decfield+decsize/2, decfield-decsize/2, decfield-decsize/2]

plt.legend(numpoints=1)

plt.savefig('M33_HSC_fields_zoom_21B.png')

#%%

m33 = sc0

fieldcen = m33

fig = plt.figure(figsize=(8,8))
ax = plt.subplot(111)
plt.ylabel(r'Dec - Dec$_0$ [deg]', fontsize=18)
and22 = SkyCoord('01h27m40.0s','28d05m25s', unit='deg')
m31 = SkyCoord('0h42m44.3s', '41d16m9s', unit='deg')

dra_and22 = (and22.ra.deg-m33.ra.value)*np.cos(m33.dec)
ddec_and22 = (and22.dec.deg-m33.dec.value)

plt.plot(dra_and22, ddec_and22, '*', color='darkviolet' ,mec='white',ms=16, label='AND XXII', zorder=5)

dp = np.loadtxt('corners_PAndAS11.txt')
ra = dp[:,0]
dec = dp[:,1]
dra = (ra-m33.ra.value)*np.cos(m33.dec)
ddec = (dec-m33.dec.value)
plt.plot(dra, ddec, lw=2, color='k')
plt.xlim(*plt.xlim()[::-1])
plt.figtext(0.43, 0.7, 'PAndAS', fontsize=20, color='k')
circle= plt.Circle((0, 0), 3.5, fill=False, color='cornflowerblue', lw=2, label='50 kpc')
circle2 = plt.Circle((0, 0), 3.5*rvir/50., fill=False, color='green', lw=2, ls=':')
circle3 = plt.Circle((0, 0), 3.5*100./50., fill=False, color='orange', lw=2, ls='--')
ax.add_artist(circle2)
ax.add_artist(circle3)
ax.add_artist(circle)
plt.xlim(12.5,-25)
plt.ylim(-12.5,25)

#######

ring = Wedge((m33.ra.deg, m33.dec.deg), 3.5*100./50., 0, 360, width=(3.5*100./50. - 3.5)),  # Full ring
p = PatchCollection(ring, alpha=0.1)
p.set_color('orange')

dra_m31 = (m31.ra.deg-m33.ra.deg)*np.cos(m33.dec)
ddec_m31 = (m31.dec.deg-m33.dec.deg)
m31ell = Ellipse((dra_m31, ddec_m31), (200./60.0), (50.0/60.0), angle=35, fill=True, linewidth=2, color='k')
ax.add_patch(m31ell)

m33ell = Ellipse((0, 0), (73./60.0), (45.0/60.0), angle=23, fill=True, linewidth=2, color='k')
ax.add_patch(m33ell)

#sc_fields_not19A_m33ref_dec = sc_fields_not19A.dec-m33.dec
#sc_fields_not19A_m33ref_ra = (sc_fields_not19A.ra-m33.ra)*np.cos(m33.dec)
#sc_fields_not19A_m33ref = SkyCoord(ra=sc_fields_not19A_m33ref_ra,dec=sc_fields_not19A_m33ref_dec,frame='icrs')

#dra_fields_not19A = (sc_fields_not19A.ra.deg-m33.ra.deg)*np.cos(m33.dec)
#ddec_fields_not19A = (sc_fields_not19A.dec.deg-m33.dec.deg)

#dra_fields19A = (sc_fields19A.ra.deg-m33.ra.deg)*np.cos(m33.dec)
#ddec_fields19A = (sc_fields19A.dec.deg-m33.dec.deg)

#m.draw_hsc_focal_planes(sc_fields_not19A_m33ref.ra.value,sc_fields_not19A_m33ref.dec.value,color='k',alpha=0.1)
#m.draw_hsc_focal_planes(sc_fields_not19A_m33ref_ra.value,sc_fields_not19A_m33ref_dec.value,color='k',alpha=0.1)

#m.draw_hsc_focal_planes(dra_fields_not19A,ddec_fields_not19A,color='k',alpha=0.1)
#m.draw_hsc_focal_planes(dra_fields19A,ddec_fields19A,color='b',alpha=0.4)

plt.legend(numpoints=1)
plt.figtext(0.57, 0.3, r'50 kpc (95-100\% coverage)', color='cornflowerblue')
plt.figtext(0.57, 0.25, r'100 kpc (49-60\% coverage)', color='orange')
plt.figtext(0.57, 0.2, r'161 kpc (35-40\% coverage)', color='green')
#plt.savefig('M33_PAndAS_sat_estimate_prelim.pdf')
plt.savefig('M33_HSC_fields_ALL_2_21B.pdf')

'''
#%%
### Output the field centers to an ASCII file:

# Create an Astropy table:
ttt = Table([sc_fields.ra.value,sc_fields.dec.value],
             names=('ra','dec'))

# write to an output file:
tmpname = 'm33_HSC_targets.txt'
ttt.write(tmpname,format='ascii',formats={'dec': '%12.6f', 'ra': '%12.6f'},overwrite=True)
print('Targets written to '+tmpname+'. Rename to prevent over-writing!')

#nn = np.repeat(targetname+'_',np.size(fieldnum_all))
#fieldnames_all0 = np.char.add(nn,np.array(fieldnum_all,dtype='str'))
#uu = np.repeat('_',np.size(fieldnum_all))
#dd = np.char.add(uu,np.array(dithnum_all,dtype='str'))
#fieldnames_all = np.char.add(fieldnames_all0,dd)
#
#ttt = Table([fieldnames_all,sc_all.ra.value,sc_all.dec.value],names=('field','ra','dec'))

# Now save only the 2020B proposed targets:
ra20B = sc_fields20B.ra.to_string(unit='hour',sep=':')
dec20B = sc_fields20B.dec.to_string(unit='degree',sep=':')
ttt20B = Table([ra20B, dec20B], names=('ra','dec'))
# write to an output file:
tmpname = 'm33_HSC_targets20B.txt'
ttt20B.write(tmpname,format='ascii',overwrite=True)
print('Targets written to '+tmpname+'. Rename to prevent over-writing!')

# Create an Astropy table:
ttt20Bd = Table([sc_fields20B.ra.value,sc_fields20B.dec.value],
             names=('ra','dec'))

# write to an output file:
tmpname = 'm33_HSC_targets20B_dec.txt'
ttt20Bd.write(tmpname,format='ascii',formats={'dec': '%12.6f', 'ra': '%12.6f'},overwrite=True)
print('Targets written to '+tmpname+'. Rename to prevent over-writing!')

'''
