# -*- coding: utf-8 -*-
"""
Created on Fri May 18 10:27:43 2018
This code compute the focal length of the OAPs and the image position.
@author: audrey.bouxin
"""
import numpy as np
import matplotlib.pyplot as plt

" Parameters"
DEG2RAD = np.pi/180;
theta_OAP0  = 20*DEG2RAD;	#[rad]  pour l'OAP0
theta_OAP1 = theta_OAP0;
theta_OAP2  = 20*DEG2RAD;	#[rad]  pour l'OAP2
#theta_OAP3  = 10*DEG2RAD;	#[rad]  pour l'OAP3
ExP2FP = 10338.74;      #[mm]
dExP   = 727.4046;      #[mm]
DDM    = 32.7;          #[mm]

# for OAP0
alpha = np.arctan(dExP*0.5/ExP2FP);  #[rad] 
PFL_OAP0 =(DDM)/2/(-1/np.tan(theta_OAP1+alpha)+(1/(np.tan(theta_OAP1+alpha))**2+1)**0.5+1/np.tan(theta_OAP1-alpha)-(1/(np.tan(theta_OAP1-alpha))**2+1)**0.5);   #[mm] pour l'OAP0
EFL_OAP0 = 2*PFL_OAP0/(1+np.cos(theta_OAP1));
print('PFL_OAP0 (estim) : ', PFL_OAP0)

 
PFL_OAP0 = 445.0;
print('PFL_OAP0 : ', PFL_OAP0)
EFL_OAP0 = 2*PFL_OAP0/(1+np.cos(theta_OAP1))
print('EFL_OAP0 : ', EFL_OAP0)
po = -(PFL_OAP0+ExP2FP);
pi = po*PFL_OAP0/(po+PFL_OAP0);
print('pi (DM)  : ', pi)

# # pour l'OAP1
fOAP1 = PFL_OAP0

 # pour l'OAP2
print('-------------------------OAP2 : ')
TTangle = 10*DEG2RAD;
DTT = 10#*np.cos(TTangle);            #[mm]
fOAP2 = DTT/DDM*fOAP1;     #[-]
alpha = np.arctan(DTT*0.5/fOAP2);    #[rad] 
PFL_OAP2 = (DTT/2)/2/(-1/np.tan(theta_OAP2+alpha)+(1/(np.tan(theta_OAP2+alpha))**2+1)**0.5+1/np.tan(theta_OAP2)-(1/(np.tan(theta_OAP2))**2+1)**0.5)  #[mm] pour l'OAP2
EFL_OAP2 = 2*PFL_OAP2/(1+np.cos(theta_OAP2))
print('PFL_OAP2 : ', PFL_OAP2)
# 
# # Pour l'OAP3
# pyr_roof = 20;  #[um]
# lambda   = 0.6; #[um]
# fnum = 60;#2*pyr_roof/lambda;
# alpha = atan(1/(fnum*2));       #[rad] pour l'OAP3
# DTT = 25.11;
# PFL_OAP3 = (DTT/2)/2/(-1/tan(theta_OAP3+alpha)+(1/(tan(theta_OAP3+alpha))**2+1)**0.5+1/tan(theta_OAP3)-(1/(tan(theta_OAP3))**2+1)**0.5)  #[mm] pour l'OAP3
# EFL_OAP3 = 2*PFL_OAP3/(1+cos(theta_OAP3))
# 




