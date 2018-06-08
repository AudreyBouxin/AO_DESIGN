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
theta_OAP0  = 15*DEG2RAD;	#[rad]  pour l'OAP0
ExP2FP = 10338.74;          #[mm]
dExP   = 727.4046;          #[mm]
DDM    = 32.7;              #[mm]

# for OAP0
alpha = np.arctan(dExP*0.5/ExP2FP);  #[rad] 
PFL_OAP0 =(DDM)/2/(-1/np.tan(theta_OAP0+alpha)+(1/(np.tan(theta_OAP0+alpha))**2+1)**0.5+1/np.tan(theta_OAP0-alpha)-(1/(np.tan(theta_OAP0-alpha))**2+1)**0.5);   #[mm] pour l'OAP0
EFL_OAP0 = 2*PFL_OAP0/(1+np.cos(theta_OAP0));
print('PFL_OAP0 (estim) : ', PFL_OAP0)

 
#PFL_OAP0 = 445.0;
print('PFL_OAP0 : ', PFL_OAP0)
EFL_OAP0 = 2*PFL_OAP0/(1+np.cos(theta_OAP0))
print('EFL_OAP0 : ', EFL_OAP0)
po = -(PFL_OAP0+ExP2FP);
pi = po*PFL_OAP0/(po+PFL_OAP0);
print('pi (DM)  : ', pi)


#Aberrations : Wabefront MAP
import zernike as zer

aa = zer.calc_zern_j(4,200,0.1,1)
print(aa)















