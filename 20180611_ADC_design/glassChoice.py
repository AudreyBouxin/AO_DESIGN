# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 09:55:36 2018

@author: audrey.bouxin
"""

from modelisationADCfunctions import glassChoiceTable
from modelisationADCfunctions import RefractiveIndex
from modelisationADCfunctions import glassPropertiesAbbeCTE3070nD
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
filenameTableResult = 'tableResults_201806201346'

if compute.lower() in ['yes','y']:
    NumberOfGlasses,tableResults = glassChoiceTable(lambda_wave)
    np.savez('results/'+filenameTableResult, NumberOfGlasses=NumberOfGlasses,tableResults=tableResults)
else:
    tmp = np.load('results/'+filenameTableResult+'.npz')
    NumberOfGlasses = tmp['NumberOfGlasses']
    tableResults = tmp['tableResults']



tableResults = tableResults.reshape(NumberOfGlasses**2,5)#tableResults.reshape(np.size(n_min)**2,5)
BDCQ_min_sorted_idx = np.argsort(tableResults[:,4])
tableResultsSorted = tableResults[BDCQ_min_sorted_idx]#glassA_idx  ,glassB_idx,  thetaA_final,  thetaB_final,  BDCQ_Final
    
    
    
plotSelection = 'plotSorted'

if plotSelection in ['plotAll','all']:
    """ Now we show the results if the keywork 'plotAll' is entered  """  
    plt.figure()
    plt.plot(tableResultsSorted[:,2]/DEG2RAD,'*',label='$\Theta_A$')
    plt.plot(tableResultsSorted[:,3]/DEG2RAD,'*',label='$\Theta_B$')
    plt.ylabel('$\Theta$ [°]')
    plt.legend()
    plt.show()
    
    plt.figure()
    plt.plot(tableResultsSorted[:,4],'*',label='BQCD')
    plt.ylabel('BQCD')
    plt.legend()
    plt.show()
    
elif plotSelection in ['plotSorted']:
#    glassA_idx  ,glassB_idx,  thetaA_final,  thetaB_final,  BDCQ_Final
    last = 1
    
#    #Plot the prism angles
#    plt.figure()    
#    plt.plot(tableResultsSorted[0:last,2]/DEG2RAD,'+',label='$\Theta_A$')
#    plt.plot(tableResultsSorted[0:last,3]/DEG2RAD,'+',label='$\Theta_B$')
#    plt.ylabel('$\Theta$ [°]')
#    plt.legend()
#    plt.show()
#    
#    #Plot the BDCQ value
#    plt.figure()
#    plt.plot(tableResultsSorted[0:last,4],'*',label='BQCD')
#    plt.ylabel('BQCD')
#    plt.legend()
#    plt.show()

    VD,CTE3070,nD = glassPropertiesAbbeCTE3070nD()
    
    """Plot the refractive index wrt lambda for the best combination of glasse"""
    lambda_um = 1e6*np.linspace(300e-9,2400e-9,100)
    n_tmp,GlassName = RefractiveIndex(lambda_um[0])

    
    for combIdx in range(last):
        print('Combination ', combIdx, ' : BDCQ = ',tableResultsSorted[combIdx,4],' | idx :', tableResultsSorted[combIdx,0], '/',tableResultsSorted[combIdx,1],
              ' | Glasses : ' , GlassName[int(tableResultsSorted[combIdx,0])],
                '/' , GlassName[int(tableResultsSorted[combIdx,1])])
        idxA = int(tableResultsSorted[combIdx,0])
        idxB = int(tableResultsSorted[combIdx,1])        
        n_A = np.array([])
        n_B = np.array([])
        for Lambda_idx in range(np.size(lambda_um)):
            n_all_glasses,GlassName = RefractiveIndex(lambda_um[Lambda_idx])
            n_A = np.append(n_A,n_all_glasses[idxA])
            n_B = np.append(n_B,n_all_glasses[idxB])


        VD_A = VD[int(tableResultsSorted[combIdx,0])]
        VD_B = VD[int(tableResultsSorted[combIdx,1])]
        nD_A = nD[int(tableResultsSorted[combIdx,0])]
        nD_B = nD[int(tableResultsSorted[combIdx,1])]
                 
        fig = plt.figure()
        ax1 = fig.add_subplot(121)
        ax1.plot(lambda_um, n_A,'+',label='n$_A$ - '+GlassName[int(tableResultsSorted[combIdx,0])])
        ax1.plot(lambda_um, n_B,'+',label='n$_B$ - '+GlassName[int(tableResultsSorted[combIdx,1])])
        ax1.set_title('Combination '+str(combIdx)+' : BDCQ = '+str(tableResultsSorted[combIdx,4])+'\n'
                  'Wedge angle : $\Theta_A$ = '+str(round(tableResultsSorted[combIdx,2]/DEG2RAD,2))+'[°] / $\Theta_B$ = '+str(round(tableResultsSorted[combIdx,3]/DEG2RAD,2))+' [°]\n'
                  'Abbe number Vd [-] : A = '+str(VD[int(tableResultsSorted[combIdx,0])])+' / B = '+str(VD[int(tableResultsSorted[combIdx,1])])+'\n'
                  'CTE-30/70 [10$^{-6}$ K$^{^-1}$]  : A = '+str(round(1e6*(CTE3070[int(tableResultsSorted[combIdx,0])]),2))+' / B = '+str(round(1e6*(CTE3070[int(tableResultsSorted[combIdx,1])]),2)))
        ax1.set_ylabel('n($\lambda$)')
        ax1.set_xlabel('$\lambda$ [um]')
        ax1.grid(True)
        ax1.legend()
        ax2 = fig.add_subplot(122)
        ax2.plot(VD_A,nD_A,'+',label=GlassName[int(tableResultsSorted[combIdx,0])])
        ax2.plot(VD_B,nD_B,'+',label=GlassName[int(tableResultsSorted[combIdx,1])])
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




