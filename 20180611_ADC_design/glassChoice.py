# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 09:55:36 2018

@author: audrey.bouxin
"""

from modelisationADCfunctions import glassChoiceTable
from modelisationADCfunctions import RefractiveIndex
#from modelisationADCfunctions import glassPropertiesAbbeCTE3070nD
import numpy as np
import matplotlib.pyplot as plt

plt.close("all")


"""
#If we want to test the computation fo the refractive index for each 
#wavelength with the corresponding Sellmeier coefficients
testRefractiveIndexCalcul()
"""

#Constants
DEG2RAD = np.pi/180.

lambda_wave = np.array([500,600,900])*1e-9  #[m]


compute = 'n'
filenameTableResult = 'tableResults_201806261438'

if compute.lower() in ['yes','y']:
    NumberOfGlasses,tableResults = glassChoiceTable(lambda_wave)
    np.savez('results/'+filenameTableResult, NumberOfGlasses=NumberOfGlasses,tableResults=tableResults)
else:
    tmp = np.load('results/'+filenameTableResult+'.npz')
    NumberOfGlasses = tmp['NumberOfGlasses']
    tableResults = tmp['tableResults']



tableResults = tableResults.reshape(NumberOfGlasses,13)
# tableResultsSorted arguments
#     0      ,      1    ,        2     ,        3     ,       4     ,         5    ,        6   ,     7  ,     8  ,       9    ,       10   ,       11     ,      12 
#glassA_idx  , glassB_idx , thetaA_final,  thetaB_final,   BDCQ_Final,        VD_A  ,       VD_B ,   nD_A ,   nD_B ,   CTE3070_A,   CTE3070_B,    GlassNameA,   GlassNameB
BDCQ_min_sorted_idx = np.argsort(tableResults[:,4].astype(float))
tableResultsSorted = tableResults[BDCQ_min_sorted_idx]


glassA_idx   = tableResultsSorted[:,0].astype(int)
glassB_idx   = tableResultsSorted[:,1].astype(int)
thetaA_final = tableResultsSorted[:,2].astype(float)
thetaB_final = tableResultsSorted[:,3].astype(float)
BDCQ_Final   = tableResultsSorted[:,4].astype(float)
VD_A         = tableResultsSorted[:,5].astype(float)
VD_B         = tableResultsSorted[:,6].astype(float)
nD_A         = tableResultsSorted[:,7].astype(float)
nD_B         = tableResultsSorted[:,8].astype(float)
CTE3070_A    = tableResultsSorted[:,9].astype(float)
CTE3070_B    = tableResultsSorted[:,10].astype(float)

    
plotSelection = 'plotSorted'

if plotSelection in ['plotAll','all']:
    """ Now we show the results if the keywork 'plotAll' is entered  """  
    #Plot the prism angles
    plt.figure(figsize=(12, 8))
    plt.plot(thetaA_final/DEG2RAD,'*',label='$\Theta_A$')
    plt.plot(thetaB_final/DEG2RAD,'*',label='$\Theta_B$')
    plt.ylabel('$\Theta$ [°]')
    plt.legend()
    plt.show()
    
    #Plot the BDCQ value
    plt.figure(figsize=(12, 8))
    plt.plot(BDCQ_Final,'*',label='BDCQ')
    plt.ylabel('BDCQ')
    plt.legend()
    plt.show()
    
elif plotSelection in ['plotSorted']:
# tableResultsSorted arguments
#     0      ,      1    ,        2     ,        3     ,       4     ,         5    ,        6   ,     7  ,     8  ,       9    ,       10   ,       11     ,      12 
#glassA_idx  , glassB_idx , thetaA_final,  thetaB_final,   BDCQ_Final,        VD_A  ,       VD_B ,   nD_A ,   nD_B ,   CTE3070_A,   CTE3070_B,    GlassNameA,   GlassNameB
    last = 10
    
#    #Plot the prism angles
#    plt.figure(figsize=(12, 8))
#    plt.plot(thetaA_final/DEG2RAD,'*',label='$\Theta_A$')
#    plt.plot(thetaB_final/DEG2RAD,'*',label='$\Theta_B$')
#    plt.ylabel('$\Theta$ [°]')
#    plt.legend()
#    plt.show()
#    
#    #Plot the BDCQ value
#    plt.figure(figsize=(12, 8))
#    plt.plot(BDCQ_Final,'*',label='BDCQ')
#    plt.ylabel('BDCQ')
#    plt.legend()
#    plt.show()


    
    """Plot the refractive index wrt lambda for the best combination of glasse"""
    lambda_um = 1e6*np.linspace(300e-9,2400e-9,100)
    

# tableResultsSorted arguments
#     0      ,      1    ,        2     ,        3     ,       4     ,         5    ,        6   ,     7  ,     8  ,       9    ,       10   ,       11     ,      12 
#glassA_idx  , glassB_idx , thetaA_final,  thetaB_final,   BDCQ_Final,        VD_A  ,       VD_B ,   nD_A ,   nD_B ,   CTE3070_A,   CTE3070_B,    GlassNameA,   GlassNameB
    
    for combIdx in range(last):
        idxA_idx       = glassA_idx[combIdx]
        idxB_idx       = glassB_idx[combIdx]
        thetaA_idx     = thetaA_final[combIdx]/DEG2RAD #[°]
        thetaB_idx     = thetaB_final[combIdx]/DEG2RAD #[°]
        BDCQ_idx       = BDCQ_Final[combIdx]
        VD_A_idx       = VD_A[combIdx]
        VD_B_idx       = VD_B[combIdx]
        nD_A_idx       = nD_A[combIdx]
        nD_B_idx       = nD_B[combIdx]
        CTE3070A_idx   = CTE3070_A[combIdx]
        CTE3070B_idx   = CTE3070_B[combIdx]
        GlassNameA_idx = tableResultsSorted[combIdx,11]
        GlassNameB_idx = tableResultsSorted[combIdx,12]
        
        print('Combination ', combIdx, ' : BDCQ = ',BDCQ_idx,' | idx :', idxA_idx, '/',idxB_idx,
              ' | Glasses : ' , GlassNameA_idx, '/' , GlassNameB_idx)
      
        n_A = np.array([])
        n_B = np.array([])
        for Lambda_idx in range(np.size(lambda_um)):
            n_all_glasses,GlassName = RefractiveIndex(lambda_um[Lambda_idx])
            n_A = np.append(n_A,n_all_glasses[idxA_idx])
            n_B = np.append(n_B,n_all_glasses[idxB_idx])
                
        fig = plt.figure(figsize=(12, 8))
        ax1 = fig.add_subplot(121)
        ax1.plot(lambda_um, n_A,'+',label='n$_A$ - '+GlassNameA_idx)
        ax1.plot(lambda_um, n_B,'+',label='n$_B$ - '+GlassNameB_idx)
        ax1.set_title('Combination '+str(combIdx)+' : BDCQ = '+str(BDCQ_idx)+'\n'
                  'Wedge angle : $\Theta_A$ = '+str(round(thetaA_idx,2))+'[°] / $\Theta_B$ = '+str(round(thetaB_idx,2))+' [°]\n'
                  'Abbe number Vd [-] : A = '+str(VD_A_idx)+' / B = '+str(VD_B_idx)+'\n'
                  'CTE-30/70 [10$^{-6}$ K$^{^-1}$]  : A = '+str(round(1e6*CTE3070A_idx,2))+' / B = '+str(round(1e6*CTE3070B_idx,2)))
        ax1.set_ylabel('n($\lambda$)')
        ax1.set_xlabel('$\lambda$ [um]')
        ax1.grid(True)
        ax1.legend()
        ax2 = fig.add_subplot(122)
        ax2.plot(VD_A_idx,nD_A_idx,'+',label=GlassNameA_idx)
        ax2.plot(VD_B_idx,nD_B_idx,'+',label=GlassNameB_idx)
        ax2.set_title('Abbe diagram')
        ax2.set_ylabel('n$_d$')
        ax2.set_xlabel('V$_D$')
        ax2.grid(True)
        ax2.legend()
        ax2.set_ylim([1.2, 2])
        ax2.set_xlim([100,20])        
        fig.show

#For example : 
# the combination cited by P.Spano during a private skype Tim8/S-FPM2 shows the following parameters :
#                     S-Tim8    S-FPM2
# Abbe Number Vd       39.24     67.73
# CTE30/70 [10^-6]      8.6      11.7
# 




