# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 08:45:41 2018

This is a function that computes the exit angle depending on the entry angle and the prism glasses

 INPUTS :
   - entryAngle   [rad] : the angle of the beam at the entrance of the prism
   - n0             [-] : The medium refractive index where the prism is (air)
   - nA,nB          [-] : Respectively the 1rst and 2nd glass refractive index
   - thetaA,thetaB [rad]: the glass angle wrt the normal of the optical axis
 OUTPUT :
   - alphaS [rad]: the exit angle

@author: jordan.voirin
"""

import numpy as np
from scipy.optimize import fmin_l_bfgs_b
import pandas as pd
import matplotlib.pyplot as plt

def exitAngle(entryAngle,n0,nA,nB,thetaA,thetaB):
    
    alpha1 = np.arcsin(n0/nA*np.sin(entryAngle))
    alpha2 = alpha1+thetaA
    alpha3 = np.arcsin(nA/nB*np.sin(alpha2))
    alpha4 = alpha3-thetaA+thetaB
    alpha5 = np.arcsin(nB/n0*np.sin(alpha4))
    alpha6 = alpha5-2*thetaB
    alpha7 = np.arcsin(n0/nB*np.sin(alpha6))
    alpha8 = alpha7+thetaB-thetaA
    alpha9 = np.arcsin(nB/nA*np.sin(alpha8))
    alpha10 = alpha9+thetaA
    alphaS = np.arcsin(nA/n0*np.sin(alpha10))
    
    return alphaS

"""
Created on Tue Jun 19 14:25:30 2018

This is a function that computes the glass refractive index

 INPUTS :
   - lambda_um [um]: the wavelength studied : HAS TO BE  a numerical value NOT ARRAY!
 OUTPUT :
   - n [-]: the refractive index array for lambda and each glass
   - glassName : the name of each glass in an array

@author: audrey.bouxin@heig-vd.ch
"""
def RefractiveIndex(lambda_um):
    SellmeierCoefficientsSchott = pd.read_csv('Glasses_catalogs/SchottSellmeierCoefficients201703.csv',header=None,names=['GlassName','B1','B2', 'B3', 'C1', 'C2', 'C3'])
#    n = (1 + SellmeierCoefficients.B1*np.power(lambda_um,2)/(np.power(lambda_um,2)-SellmeierCoefficients.C1) 
#            + SellmeierCoefficients.B2*np.power(lambda_um,2)/(np.power(lambda_um,2)-SellmeierCoefficients.C2) 
#            + SellmeierCoefficients.B3*np.power(lambda_um,2)/(np.power(lambda_um,2)-SellmeierCoefficients.C3))**0.5;
    n_schott = ((1 + np.array(SellmeierCoefficientsSchott.B1)*np.power(lambda_um,2)/(np.power(lambda_um,2)-np.array(SellmeierCoefficientsSchott.C1)) 
        + np.array(SellmeierCoefficientsSchott.B2)*np.power(lambda_um,2)/(np.power(lambda_um,2)-np.array(SellmeierCoefficientsSchott.C2)) 
        + np.array(SellmeierCoefficientsSchott.B3)*np.power(lambda_um,2)/(np.power(lambda_um,2)-np.array(SellmeierCoefficientsSchott.C3))))**0.5
    
    SellmeierCoefficientsOhara = pd.read_csv('Glasses_catalogs/OharaSellmeierCoefficients201805.csv',header=None,names=['GlassName','B1','B2', 'B3', 'C1', 'C2', 'C3'])
    n_ohara = ((1 + np.array(SellmeierCoefficientsOhara.B1)*np.power(lambda_um,2)/(np.power(lambda_um,2)-np.array(SellmeierCoefficientsOhara.C1)) 
        + np.array(SellmeierCoefficientsOhara.B2)*np.power(lambda_um,2)/(np.power(lambda_um,2)-np.array(SellmeierCoefficientsOhara.C2)) 
        + np.array(SellmeierCoefficientsOhara.B3)*np.power(lambda_um,2)/(np.power(lambda_um,2)-np.array(SellmeierCoefficientsOhara.C3))))**0.5
    
    n = np.concatenate((n_schott, n_ohara),axis=0)

    glassName = np.concatenate((SellmeierCoefficientsSchott.GlassName, SellmeierCoefficientsOhara.GlassName),axis=0)
    return n,glassName

"""
Created on Tue Jun 20 12:46:32 2018

This is a function tests the refractive index calculs with the Sellmeier coefficients

 INPUTS : [-]
 OUTPUT : [-]

@author: audrey.bouxin@heig-vd.ch
"""
def testRefractiveIndexCalcul():
    lambda_wave = 1e-9*np.linspace(300,2500,100)
    lambda_wave_um = lambda_wave*1e6
    n = np.array([])
    glassName =  np.array([])
    for idx in range(np.size(lambda_wave)):
        n_tmp,glassName = RefractiveIndex(lambda_wave_um[idx])
        n = np.append(n,n_tmp)
    
    n = n.reshape(np.size(lambda_wave),5)
    plt.figure()
    plt.hold(True)
    for idx in range(5):
        plt.plot(lambda_wave_um*1e3,n[:,idx],label=glassName[idx])
    plt.xlabel('\lambda [nm]')
    plt.ylabel('n')
    plt.legend()
    axes = plt.gca()
    axes.set_xlim([lambda_wave[0]*1e9, lambda_wave[-1]*1e9])
    axes.set_ylim([1.1,1.9])
    plt.show()
    
    return 0

# -*- coding: utf-8 -*-
"""
This is a function that computes matrix containing all the glasses combinations performance.

 INPUTS :
   - lambda_wave  [m] Studied wavelengths (can be an array)
 OUTPUT :
   - tableResults  the matrix containing the glass combination performance, 
                   the angles giving the best geometrical results for each prims pair and the metric value
   
Created on Tue Jun 19 09:16:04 2018

@author: audrey.bouxin (audrey.bouxin@heig-vd.ch)
 
 June 2018, ABx Creation 
"""    
    
def glassChoiceTable(lambda_wave):
    #Constants
    DEG2RAD = np.pi/180.
    lambda_um = lambda_wave*1e6            #[um]

    n_min,glassName_min = RefractiveIndex(lambda_um[0])
    n_cen,glassName_cen = RefractiveIndex(lambda_um[1])
    n_max,glassName_max = RefractiveIndex(lambda_um[2])
    GlassName = glassName_min
    
    Ratm = Refraction_atmosphere(lambda_wave, 70,standard=0)
    entryAngles    = Ratm-Ratm[1]
    entryAngle_min = entryAngles[0]
    entryAngle_av  = entryAngles[1]
    entryAngle_max = entryAngles[2]
    n0 = 1.
    thetaMin = -20*DEG2RAD
    thetaMax = 20*DEG2RAD
    
    tableResults = np.array([])
    
    def func2min_BDCQ(X,*args):
        thetaA=X[0]
        thetaB=X[1]
        exitAngle_min = exitAngle(entryAngle_min,n0,n_min[glassA_idx],n_min[glassB_idx],thetaA,thetaB)
        exitAngle_av  = exitAngle(entryAngle_av ,n0,n_cen[glassA_idx],n_cen[glassB_idx],thetaA,thetaB)
        exitAngle_max = exitAngle(entryAngle_max,n0,n_max[glassA_idx],n_max[glassB_idx],thetaA,thetaB)
        BDCQ = ((exitAngle_min-exitAngle_av)**2 + (exitAngle_max-exitAngle_av)**2 + (exitAngle_av-entryAngle_av)**2)**.5 #Beam Dispersion Correction Quality metric
        return BDCQ
    
    NumberOfGlasses=np.size(n_min)
        
    
    for glassA_idx in range(NumberOfGlasses):
        for glassB_idx in range(NumberOfGlasses):
            print('Num glass A : ',glassA_idx, ' Glass IDs : ', ' A : ',GlassName[glassA_idx])#,' B : ',GlassName[glassB_idx])
            initial_values = np.array([thetaMin, thetaMin])    #glassA_idx, glassB_idx, thetaA, thetaB
            bounds = [(thetaMin,thetaMax),(thetaMin,thetaMax)]
            thetas_final,BDCQ_Final,info = fmin_l_bfgs_b(func=func2min_BDCQ, x0=initial_values, fprime=None, args=(glassA_idx,glassB_idx), approx_grad=True, bounds=bounds, m=10, factr=10000000.0, pgtol=1e-05, epsilon=1e-08, iprint=-1, maxfun=15000, maxiter=15000, disp=None, callback=None, maxls=20)
            thetaA_final = thetas_final[0]#[rad]
            thetaB_final = thetas_final[1]#[rad]
            tableResults = np.append(tableResults,np.array([glassA_idx,glassB_idx,thetaA_final,thetaB_final,BDCQ_Final]))
            
    return NumberOfGlasses,tableResults
   
# -*- coding: utf-8 -*-
"""
This is a function that looks for the Abbe number in the glass properties file

 INPUTS :
   - lambda_wave  [m] Studied wavelengths (can be an array)
   - zenith_angle [°] Zenith angle of observation
 OUTPUT :
   - Ratm       [rad] angular dispersion of the atmosphere
   
Created on Tue Jun 12 09:16:04 2018

@author: audrey.bouxin (audrey.bouxin@heig-vd.ch)
 
 June 2018, ABx Creation 
"""

def glassPropertiesAbbeCTE3070nD():
    #glassA_idx  ,glassB_idx,  thetaA_final,  thetaB_final,  BDCQ_Final
    schottGlassProperties = pd.read_csv('Glasses_catalogs/SchottGlassProperties201703.csv',header=None,names=['VD','alpha3070','alpha20300','nD'])
    oharaGlassProperties = pd.read_csv('Glasses_catalogs/OharaGlassProperties201703.csv',header=None,names=['VD','alpha3070','alpha100300','nD'])
    
    VD = np.concatenate((schottGlassProperties.VD, oharaGlassProperties.VD),axis=0)
    CTE3070 = np.concatenate((schottGlassProperties.alpha3070*1e-6, oharaGlassProperties.alpha3070*1e-7),axis=0)
    nD = np.concatenate((schottGlassProperties.nD, oharaGlassProperties.nD),axis=0)
    return VD,CTE3070,nD





# -*- coding: utf-8 -*-
"""
This is a function that computes the atmosphere refraction angle at specified wavelengths.

 INPUTS :
   - lambda_wave  [m] Studied wavelengths (can be an array)
   - zenith_angle [°] Zenith angle of observation
 OUTPUT :
   - Ratm       [rad] angular dispersion of the atmosphere
   
Created on Tue Jun 12 09:16:04 2018

@author: audrey.bouxin (audrey.bouxin@heig-vd.ch)
 
 June 2018, ABx Creation 
"""

def Refraction_atmosphere(lambda_wave, zenith_angle, standard):

    #constants:
    DEG2RAD = np.pi/180.
    zenith_angle_rad = DEG2RAD*zenith_angle;
    lambda_wave_um = lambda_wave*1e6;
    g  = 9.80665;        #[m/s2]  earth-surface gravitational acceleration
    Mw = 18.01528e-3;    #[kg/mol] (Ciddor 1996)
    R  = 8.314510;       #[J/mol/K] (Davis 1992 & Ciddor 1996)
    
    if standard ==1:   
        #General parameters
        Height = 0.             #[m] altitude standard
        T_C    = 15.    ;       #[°C] temperature
        T_K    = T_C + 273.15;  #[K] temperature
        T0_C   = 15.            #[°C] standard temperature (Robo-AO)
        T0_K   = T0_C+273.15;   #[K] standard temperature 
        DT     = 0.0065;        #[K] vertical gradient of temperature 0.65K for 100m (Wikipedia)
        P0     = 1.01325e5;     #[Pa] Normal pressure at altitude = 0
        RH     = 0.5            #[-] relative humidity       
    else:
                #General parameters
        Height = 3170;          #[m] altitude Karakaya
        T_C    = -15;           #[°C] temperature
        T_K    = T_C + 273.15;  #[K] temperature
        T0_C   = 15.            #[°C] standard temperature (Robo-AO)
        T0_K   = T0_C+273.15;   #[K] standard temperature 
        DT     = 0.0065;        #[K] vertical gradient of temperature 0.65K for 100m (Wikipedia)
        P0     = 1.01325e5;     #[Pa] Normal pressure at altitude = 0
        RH     = 0.8;           #[-] relative humidity     
    
    
    
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
    #[rad] angular dispersion of the atmosphere
    Ratm = kappa*gamma*(1-beta)*np.tan(zenith_angle_rad)-kappa*gamma*(beta-gamma/2)*(np.tan(zenith_angle_rad))**3;    
    
    return Ratm

























