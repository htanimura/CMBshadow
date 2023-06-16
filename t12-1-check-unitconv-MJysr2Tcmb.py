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

# Planck frequencies
freq = np.array([100, 143, 217, 353, 545, 857])
Nf = len(freq)

for I in range(Nf):
    freq_i = freq[I] * u.GHz
    equiv = u.thermodynamic_temperature(freq_i, Planck18.Tcmb0)

    print(freq_i, (1. * u.K).to(u.MJy / u.sr, equivalencies=equiv))

#U_C [MJy/sr/KCMB]
#100-avg	244.0960 ± 0.3050
#143-avg	371.7327 ± 0.0784	
#217-avg	483.6874 ± 0.0118
#353-avg	287.4517 ± 0.0085
#545-avg	58.0356 ± 0.0278
#857-avg	2.2681 ± 0.0270


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


