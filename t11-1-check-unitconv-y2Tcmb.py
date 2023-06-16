#! /usr/bin/python
# coding: utf-8

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from astropy.constants import c, m_e, k_B, h, sigma_T, kpc, pc
from astropy.constants import codata2018 as const 
#print(const.h)
import astropy.units as u
from astropy.cosmology import Planck18, z_at_value
#z_at_value(Planck18.age, 2 * u.Gyr)

Nside = 512
indir = "./data/"
simdir = "./data/observations/litebird/"
Tcmb0 = Planck18.Tcmb0.value #2.7255 #[K]

freq = np.array([40,50,60,68,78,89,100,119,140,166,195,235,280,337,402]) #GHz
freq_char = ['040','050','060','068','078','089','100','119','140','166','195','235','280','337','402'] #GHz
beam_fwhm = np.array([70.5, 58.5, 51.1, 41.6, 36.9, 33.0, 30.2, 26.3, 23.7, 28.9, 28.0, 24.7, 22.5, 20.9, 17.9]) # arcmin in FWHM
beam_sigma = beam_fwhm/sqrt(8*log(2))

# Read data

#for ID in range(15):
for ID in [4]:
    #ID = 4
    freq_i = freq[ID]
    fwhm_i = beam_fwhm[ID]; sigma_i = beam_sigma[ID]
    print("freq: %d GHz, beam (FWHM): %.1f arcmin" %(freq_i, fwhm_i))
    
    ymap = hp.read_map(indir+'thermalsz_map.fits')
    fmap = hp.read_map(simdir+'%dGHz/group3_map_%dGHz.fits' %(freq_i, freq_i))
    
    v = np.array([freq_i]) * 1e9 # GHz -> Hz
    x = (h.value*v)/(k_B.value*Tcmb0)
    
    fx = x/np.tanh(x/2.) - 4
    
    dT_tsz = Tcmb0 * ymap * fx * 1e6 # y -> uK_Tcmb
    dT_tsz_beam = hp.smoothing(dT_tsz, sigma=np.deg2rad(sigma_i/60.))

    """
    plt.figure()
    plt.plot(fmap, dT_tsz_beam, ',')
    plt.xlabel('$T_{%d GHz}$ [uK_CMB] from LiteBIRD sim' %(freq_i))
    plt.ylabel('$T_{%d GHz}$ [uK_CMB] from thermalsz_map.fits' %(freq_i))
    plt.xlim(-50,0)
    plt.ylim(-50,0)
    """    
    plt.figure()
    plt.hist(dT_tsz_beam/fmap, bins=100)
    plt.xlabel('My/LiteBIRD sim (%d GHz)' %(freq_i))




"""
### thermodynamic temperature
from astropy import units as u
from astropy.cosmology import Planck15
freq = 143 * u.GHz
equiv = u.thermodynamic_temperature(freq, Planck15.Tcmb0)
print((1. * u.mK).to(u.MJy / u.sr, equivalencies=equiv)  )

### brightness temperature
beam_sigma = 50*u.arcsec
beam_area = 2*np.pi*(beam_sigma)**2
freq = 5*u.GHz
equiv = u.brightness_temperature(freq)
print((1*u.Jy/beam_area).to(u.K, equivalencies=equiv))
"""


