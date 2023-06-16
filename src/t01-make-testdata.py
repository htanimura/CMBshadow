import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from cosmology import *
import astropy.units as u
from astropy.cosmology import Planck18
from astropy.constants import c, m_e, k_B, h, sigma_T, kpc, pc
from astropy.constants import codata2018 as const 
from astropy.io import fits
from astropy.coordinates import SkyCoord
from hpproj import hp_project

simdir = '/gpfs02/work/tanimura/data/sims/websky/'

Nside = 2048; Npix = hp.nside2npix(Nside);
Tcmb0 = Planck18.Tcmb0.value #2.7255 #[K]

### Planck configuration
freq = np.array([100, 143, 217, 353, 545, 857]) # GHz
Nf = len(freq)
beam = np.array([9.68, 7.30, 5.02, 4.94, 4.83, 4.64]) # [arcmin] in FWHM
noise = np.array([32.27738969, 13.77801944, 19.36806465, 60.44045947, 342.99701177, 8442.34293558]) #[uK]

### read data

# Unit conversion MJy/sr -> Tcmb [K]
Kcmb2MJy = np.zeros(Nf)
for I in range(Nf):
    freq_i = freq[I] * u.GHz
    equiv = u.thermodynamic_temperature(freq_i, Planck18.Tcmb0)
    #print(freq_i, (1. * u.K).to(u.MJy / u.sr, equivalencies=equiv))
    Kcmb2MJy[I] = (1. * u.K).to(u.MJy / u.sr, equivalencies=equiv).value


# CMB in uK
cmb_alm = hp.read_alm(simdir + 'lensed_alm.fits') # [uK]
cmb_map = hp.alm2map(cmb_alm, Nside)

#cmb_cl = hp.alm2cl(cmb_alm)
#ell = np.arange(len(cmb_cl))
#plt.plot(ell, ell*(ell+1)/(2.*np.pi) * cmb_cl)
#plt.xlim(0,2500)

# CIB in uK
cibs = np.zeros((Nf, Npix))

for I in range(Nf):
    cib_fi_4096 = hp.read_map(simdir + 'cib_nu0100.fits') # [MJy/sr]
    cib_fi_2048 = hp.ud_grade(cib_fi_4096, Nside)
    cibs[I] = cib_fi_2048 / Kcmb2MJy[I] * 1e6 # [MJy/sr] -> [uK]

# tSZ in uK
tsz = hp.read_map(simdir + 'tsz_2048.fits')

v = freq * 1e9 # GHz -> Hz
x = (h.value*v)/(k_B.value*Tcmb0)
fx = x/np.tanh(x/2.) - 4

tszs = np.zeros((Nf, Npix))
for I in range(Nf):
    tszs[I] = Tcmb0 * tsz * fx[I] * 1e6 # y -> [uK]

# noise in uK
np.random.seed(seed=32)
noises = np.zeros((Nf, Npix))
for I in range(Nf):
    noises[I] = numpy.random.normal(loc=0, scale=noise[I], size=Npix)




# Halo catalog
rho = 2.775e11*omegam*h**2 # Msun/Mpc^3
f=open(simdir + 'halos.pksc')
N=np.fromfile(f,count=3,dtype=np.int32)[0]

# only take first N entries for testing (there are ~8e8 halos total...)
# comment the following line to read in all halos
N = 3

catalog=np.fromfile(f,count=N*10,dtype=np.float32)
catalog=np.reshape(catalog,(N,10))

x  = catalog[:,0];  y = catalog[:,1];  z = catalog[:,2] # Mpc (comoving)
vx = catalog[:,3]; vy = catalog[:,4]; vz = catalog[:,5] # km/sec
R  = catalog[:,6] # Mpc

# convert to mass, comoving distance, radial velocity, redshfit, RA and DEc
M200m    = 4*np.pi/3.*rho*R**3        # this is M200m (mean density 200 times mean) in Msun
chi      = np.sqrt(x**2+y**2+z**2)    # Mpc
vrad     = (x*vx + y*vy + z*vz) / chi # km/sec
redshift = zofchi(chi)      

#theta, phi  = hp.vec2ang(np.column_stack((x,y,z))) # in radians
glon, glat  = hp.vec2ang(np.column_stack((x,y,z)), lonlat=True) # in degrees



# Plot (healpix to 2d projection)
#hp_data, hp_header = hp.read_map(simdir+'tsz_2048.fits', h=True, nest=None)
hp_header = {'NSIDE': Nside,
'ORDERING': 'RING',
'COORDSYS': 'G'}
hp_hdu = fits.ImageHDU(tsz, fits.Header(hp_header))

J = 0
pixsize_deg = 1./60. 
coord = SkyCoord(glon[J], glat[J], unit="deg", frame="galactic")
hdu = hp_project(hp_hdu, coord, pixsize=pixsize_deg, shape_out=(60, 60))
plt.imshow(hdu.data, origin='lower', interpolation='none')


