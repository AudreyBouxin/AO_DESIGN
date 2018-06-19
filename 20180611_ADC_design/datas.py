# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 11:38:15 2018

@author: audrey.bouxin
"""

import pandas as pd
from modelisationADCfunctions import Refraction_atmosphere
from modelisationADCfunctions import exitAngle
import numpy as np
import matplotlib.pyplot as plt

#Constants
DEG2RAD = np.pi/180.

lambda_wave = np.array([500,600,900])*1e-9  #[m]
lambda_wave_um = lambda_wave*1e6            #[um]
SellmeierCoefficients = pd.read_csv('Glasses_catalogs/SellmeierCoefficients.csv',header=None,names=['GlassName','B1','C1','B2','C2','B3','C3'])

lambda_min = lambda_wave_um[0]
n_min = (1 + SellmeierCoefficients.B1*lambda_min**2/(lambda_min**2-SellmeierCoefficients.C1) 
    + SellmeierCoefficients.B2*lambda_min**2/(lambda_min**2-SellmeierCoefficients.C2) 
    + SellmeierCoefficients.B3*lambda_min**2/(lambda_min**2-SellmeierCoefficients.C3))**0.5;

lambda_av = lambda_wave_um[1]
n_av = (1 + SellmeierCoefficients.B1*lambda_av**2/(lambda_av**2-SellmeierCoefficients.C1) 
    + SellmeierCoefficients.B2*lambda_av**2/(lambda_av**2-SellmeierCoefficients.C2) 
    + SellmeierCoefficients.B3*lambda_av**2/(lambda_av**2-SellmeierCoefficients.C3))**0.5;

lambda_max = lambda_wave_um[2]
n_max = (1 + SellmeierCoefficients.B1*lambda_max**2/(lambda_max**2-SellmeierCoefficients.C1) 
    + SellmeierCoefficients.B2*lambda_max**2/(lambda_max**2-SellmeierCoefficients.C2) 
    + SellmeierCoefficients.B3*lambda_max**2/(lambda_max**2-SellmeierCoefficients.C3))**0.5;

thetaMin = -10
thetaMax = 10
Nbpt = (thetaMax-thetaMin)/0.1+1
thetaA = np.linspace(thetaMin, thetaMax,Nbpt)*DEG2RAD
thetaB = thetaA


Ratm = Refraction_atmosphere(lambda_wave, 70)
entryAngles = Ratm-Ratm[1]
entryAngle_min = entryAngles[0]
entryAngle_av  = entryAngles[1]
entryAngle_max = entryAngles[2]
n0 = 1.

tableResults = np.array([])

for glassA_idx in range(np.size(n_min)):
    for glassB_idx in range(np.size(n_min)):
        if glassA_idx ==glassB_idx: continue
        for thetaA_idx in range(int(Nbpt)):
            for thetaB_idx in range(int(Nbpt)):
                exitAngle_min = exitAngle(entryAngle_min,n0,n_min[glassA_idx],n_min[glassB_idx],thetaA[thetaA_idx],thetaB[thetaB_idx])
                exitAngle_av  = exitAngle(entryAngle_av ,n0,n_av[glassA_idx],n_av[glassB_idx], thetaA[thetaA_idx],thetaB[thetaB_idx])
                exitAngle_max = exitAngle(entryAngle_max,n0,n_max[glassA_idx],n_max[glassB_idx],thetaA[thetaA_idx],thetaB[thetaB_idx])
                BDCQ = (exitAngle_min-exitAngle_av)**2 + (exitAngle_max-exitAngle_av)**2 + (exitAngle_av-entryAngle_av)**2 #Beam Dispersion Correction Quality metric
                tableResults = np.append(tableResults,np.array([glassA_idx,glassB_idx, thetaA[thetaA_idx],thetaB[thetaB_idx],BDCQ]))


BDCQ_min_sorted_idx = np.argsort(tableResults[:,4])








