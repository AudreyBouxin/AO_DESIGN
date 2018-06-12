# -*- coding: utf-8 -*-
"""
Created on Mon Jun 11 11:36:34 2018

Code to test the refraction_atmosphere function

@author: audrey.bouxin
"""
from Refraction_atmosphere import Refraction_atmosphere
import numpy as np
import matplotlib.pyplot as plt

RAD2ASEC = 3600*180/np.pi

# =============================================================================
#  INPUTS :
#    - lambda_wave  [m] Studied wavelengths
#    - zenith_angle [°] Zenith angle of observation
#  OUTPUT :
#    - Ratm       [rad] la dispersion angulaire de l'atmosphère
# =============================================================================
lambda_wave = np.linspace(500,900,num=100)*1e-9
zenith_angle = 45
test_Ratm = Refraction_atmosphere(lambda_wave,zenith_angle)
plt.figure()
plt.plot(lambda_wave*1e9, test_Ratm*RAD2ASEC)
plt.ylabel('Refraction R($\lambda$) [asec]')
plt.xlabel('Wavelength [nm]')
plt.axis([lambda_wave[0]*1e9, lambda_wave[-1]*1e9, min(test_Ratm*RAD2ASEC),max(test_Ratm*RAD2ASEC)])
plt.show()