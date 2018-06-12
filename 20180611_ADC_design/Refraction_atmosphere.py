# -*- coding: utf-8 -*-
"""
This is a function that computes the atmosphere refraction angle at specified wavelengths.

 INPUTS :
   - lambda_wave  [m] Studied wavelengths
   - zenith_angle [°] Zenith angle of observation
 OUTPUT :
   - Ratm       [rad] la dispersion angulaire de l'atmosphère
   
Created on Tue Jun 12 09:16:04 2018

@author: audrey.bouxin (audrey.bouxin@heig-vd.ch)
 
 June 2018, ABx Creation 
"""

import numpy as np


def Refraction_atmosphere(lambda_wave, zenith_angle):

    #constants:
    DEG2RAD = np.pi/180.
    zenith_angle_rad = DEG2RAD*zenith_angle;
    lambda_wave_um = lambda_wave*1e6;
    g  = 9.80665;        #[m/s2]  earth-surface gravitational acceleration
    Mw = 18.01528e-3;    #[kg/mol] (Ciddor 1996)
    R  = 8.314510;       #[J/mol/K] (Davis 1992 & Ciddor 1996)
    
    #General parameters
    Height = 0.#3170;       #[m] altitude Karakaya
    T_C    = 15.#-15;       #[°C] temperature
    T_K    = T_C + 273.15;  #[K] temperature
    T0_C   = 15.            #[°C] standard temperature (Robo-AO)
    T0_K   = T0_C+273.15;   #[K] standard temperature 
    DT     = 0.0065;        #[K] vertical gradient of temperature 0.65K for 100m (Wikipedia)
    P0     = 1.01325e5;     #[Pa] Normal pressure at altitude = 0
    RH     = 0.5;#♣0.8;     #[-] relative humidity       
    
    
    
    
# =============================================================================
#     Point 1 Ciddor : svp, f, xw
# =============================================================================
    
    ####### Saturated water vapor pressure #######
    #Coefficients for the saturated vapor pressure svp (Davis 1992 & Ciddor 1996)
    A = 1.2378847e-5;   #[K^-2]
    B = -1.9121316e-2;  #[K^-1]
    C = 33.93711047;    #[-]
    D = -6.4341645e3;   #[K]
    #pressure of the saturated water vapor in the humid air (Davis 1992 & Ciddor 1996)
    svp = np.exp(A*T_K**2+B*T_K+C+D/T_K);    #[Pa]       
    
    
    ####### Enhancement factor of water vapor in the air #######
    # =============================================================================
    #     Point 3 Ciddor : Ma
    # ============================================================================= 
    xCO2 = 400;          #[ppm] (Davis 1992 & Ciddor 1996)
    Ma = (28.9635+12.011*1e-6*(xCO2-400))*1e-3;   #[kg/mol] molar mass of dry air containing xCO2 ppm of CO2 (Davis 1992 & Ciddor 1996)
    P  = P0*(1-DT*Height/T0_K)**(g*Ma/(R*DT));  #[Pa] the total pressure (Wikipedia)
    #Coefficients for f (Davis 1992 & Ciddor 1996)
    alpha = 1.00062;    #[-]
    beta  = 3.14e-8;    #[Pa^-1]
    gamma = 5.6e-7;     #[°C^-2]
    #Enhancement factor of water vapor in the air (Davis 1992 & Ciddor 1996)
    f = alpha+beta*P+gamma*T_C**2;  #[-]
    
    
    ####### Molar fraction of water vapor in moist air #######
    xw = f*RH*svp/P  #[-] (Ciddor 1996) Note : the partial vapor pw = RH*svp
    
    
# =============================================================================
#     Point 2 Ciddor : n_axs
# =============================================================================   
    #Constants involved in the standard phase and group refractivities of dry air (Ciddor 1996)
    k0 = 238.0185;  #[um^-2]
    k1 = 5792105;   #[um^-2]
    k2 = 57.362;    #[um^-2]
    k3 = 167917;    #[um^-2]
    
    #Equations for the atmospheric refraction index in the visible and near infrared
    sigma = 1/lambda_wave_um;  #[um^-1] wave number

    #refractive index of standard air at 15°C, 101325Pa,0#humidity, 450ppm of CO2
    n_as = 1+1e-8*(k1/(k0-sigma**2)+k3/(k2-sigma**2));  
    
    #For us it is not usefull to have the following equation because we are in
    #a standard case so naxs=nas, but I let it in order to remember the sequence
    #of equations
    xc = 450.;                                   #[ppm] number of ppm of CO2 (standard) (Ciddar 1996)
    n_axs = 1+(n_as-1)*(1+0.534e-6*(xc-450));    #refractive index of standard air at 15°C, 101325Pa,0%humidity, xc ppm of CO2
    

# =============================================================================
#     Point 4 Ciddor : Za
# =============================================================================  
    #Coefficients the compressibility factor Z (Davis 1992 & Ciddor 1996)
    a0 = 1.58123e-6; #[K*Pa^-1]
    a1 = -2.9331e-8; #[Pa^-1]
    a2 = 1.1043e-10; #[(K*Pa)^-1]
    b0 = 5.707e-6;   #[K*Pa^-1]
    b1 = -2.051e-8;  #[Pa^-1]
    c0 = 1.9898e-4;  #[K*Pa^-1]
    c1 = -2.376e-6;  #[Pa^-1]
    d  = 1.83e-11;   #[K^2Pa^-2]
    e  = -0.765e-8;  #[K^2Pa^-2]     
    
    ####### The compressibility of dry air #######
    xw_as = 0  #[-] Molar fraction of water vapor in moist air : here we are in dry air so 0
    P_as = P0       #[Pa] (Ciddor 1996)
    T_as = T0_K     #[K] (Ciddor 1996)
    Zas = 1-P_as/T_as*(a0+a1*T_as+a2*T_as**2+(b0+b1*T_as)*xw_as +(c0+c1*T_as)*xw_as**2) +(d+e*xw_as**2)*(P0/T_as)**2;   #[-] compressibility of dry air (Davis 1992 & Ciddor 1996)
    
# =============================================================================
#     Point 5 Ciddor : Zw
# =============================================================================    
    ####### The compressibility of pure water vapor #######
    xw_ws = 1  #[-] Molar fraction of water vapor in moist air : here we are pure water vapor so 1
    P_ws = 1333         #[Pa] (Ciddor 1996)
    T_ws = 273.15+20    #[K] (Ciddor 1996)
    Zws = 1-P_ws/T_ws*(a0+a1*T_ws+a2*T_ws**2+(b0+b1*T_ws)*xw_ws +(c0+c1*T_ws)*xw_ws**2) +(d+e*xw_ws**2)*(P_ws/T_ws)**2;   #[-] compressibility of dry air (Davis 1992 & Ciddor 1996)
    
    
# =============================================================================
#     Point 6 Ciddor : ro_axs, ro_ws
# =============================================================================     
    ####### The density of standard air #######
    ro_axs = (P_as*Ma/(Zas*R*T_as))*(1-xw_as*(1-Mw/Ma))  #[kg/m^3] The density of standard air (Ciddor 1996)
    
    ####### The density of standard air #######
    ro_ws = (P_ws*Ma/(Zws*R*T_ws))*(1-xw_ws*(1-Mw/Ma))  #[kg/m^3] The density of standard water vapor (Ciddor 1996)
    
    
# =============================================================================
#     Point 7 Ciddor : Z
# =============================================================================    
    ####### The compressibility of moist air under experimental conditions #######
    xw_w = xw   #[-] cf Point 1 (Ciddor 1996)
    P_w  = P    #[Pa] 
    T_w  = T_K  #[K] 
    Z   = 1-P_w/T_w*(a0+a1*T_w+a2*T_w**2+(b0+b1*T_w)*xw_w +(c0+c1*T_w)*xw_w**2) +(d+e*xw_w**2)*(P_w/T_w)**2;   #[-] compressibility of moist air (Davis 1992 & Ciddor 1996)    
    
# =============================================================================
#     Point 8 Ciddor : ro_a
# =============================================================================      
    ####### The density of the dry component of the moist air #######
    P_a = P   #[Pa] the total atmospheric pressure
    xa = xw   #[-] cf Point 1 (Ciddor 1996)
    T_a = T_K #[K] (Ciddor 1996)
    ro_a = P_a*Ma*(1-xa)/(Z*R*T_a)  #[kg/m^3] The density of the dry component of the moist air (Ciddor 1996)
    
# =============================================================================
#     Point 9 Ciddor : ro_w
# =============================================================================      
    ####### The density of the dry component of the moist air #######
    P_w = P   #[Pa] the total atmospheric pressure
    xw = xw   #[-] cf Point 1 (Ciddor 1996)
    T_w = T_K #[K] (Ciddor 1996)
    ro_w = P_w*Mw*xw/(Z*R*T_w)  #[kg/m^3] The density of the dry component of the moist air (Ciddor 1996)    
    
# =============================================================================
#     Point 10 Ciddor : n_prop
# =============================================================================   
    ####### The refractive index of water vapor at standard conditions #######
    #Constants involved in the standard phase and group refractivities of water vapor (Ciddor 1996)
    w0 = 295.235;   #[um^-2]
    w1 = 2.6422;    #[um^-2]
    w2 = -0.032380;	#[um^-4]
    w3 = 0.004028;  #[um^-6]
    #refractive index of water vapor at 20°C, 1335Pa
    n_ws = 1+1.022e-8*(w0+w1*sigma**2 +w2*sigma**4+w3*sigma**6); #refractive index of water vapor at 20°C, 1335Pa (Ciddor 1996)


    ####### The refractive index of moist air (atmosphere) #######
    n_prop = 1+(ro_a/ro_axs)*(n_axs-1)+(ro_w/ro_ws)*(n_ws-1); #(Ciddor 1996)   
    
    
###############################################################################    
############################################################################### 
############################################################################### 
    
    
# =============================================================================
#    Calcul of the atmospheric dispersion angle
# =============================================================================  
    kappa = 1.;  #constant for a spherical Earth kappa = 1 (Stone 1996)
    beta = 0.001254*(T_K/273.15);   #(Stone 1996)
    n = n_prop
    gamma = (n-1)   # #(Stone 1996)
    #we can retain only the first two terms of the expanded power serie of R for zenith angles < 75° (Stone 1996)
    Ratm = kappa*gamma*(1-beta)*np.tan(zenith_angle_rad)-kappa*gamma*(beta-gamma/2)*(np.tan(zenith_angle_rad))**3;    
    
    return Ratm
    
    
    
    