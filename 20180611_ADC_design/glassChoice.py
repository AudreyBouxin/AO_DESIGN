# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 09:55:36 2018

@author: audrey.bouxin
"""

import pandas as pd
#from modelisationADCfunctions import Refraction_atmosphere
#from modelisationADCfunctions import exitAngle
from modelisationADCfunctions import sellmeierRefractiveIndex
from modelisationADCfunctions import glassChoiceTable
import numpy as np
import matplotlib.pyplot as plt
plt.close("all")


#Constants
DEG2RAD = np.pi/180.

lambda_wave = np.array([500,600,900])*1e-9  #[m]

SellmeierCoefficients_Schott = pd.read_csv('Glasses_catalogs/SellmeierCoeffSchott.csv',header=None,names=['GlassName','B1','C1','B2','C2','B3','C3'])
DispersionCoefficients_Ohara = pd.read_csv('Glasses_catalogs/SellmeierCoeffOhara.csv', header=None,names=['GlassName','A1','A2','A3','A4','A5','A6'])

#NumberOfGlasses,tableResults = glassChoiceTable(lambda_wave,SellmeierCoefficients)
#np.save('tableResults201806191417', NumberOfGlasses,tableResults)

NumberOfGlasses = np.size(SellmeierCoefficients.B1)
tableResults = np.load('tableResults201806191013.npy')
tableResults = tableResults.reshape(NumberOfGlasses**2,5)#tableResults.reshape(np.size(n_min)**2,5)
BDCQ_min_sorted_idx = np.argsort(tableResults[:,4])
tableResultsSorted = tableResults[BDCQ_min_sorted_idx]
    
    
    
plotSelection = 'plotSorted10'

if plotSelection in ['ploAll','all']:
    """ Now we show the results if the keywork 'plotAll' is entered  """  
    plt.figure()
    plt.plot(tableResultsSorted[:,2]/DEG2RAD,'*',label='theta A')
    plt.plot(tableResultsSorted[:,3]/DEG2RAD,'*',label='theta B')
    plt.ylabel('theta [°]')
    plt.legend()
    plt.show()
    
    plt.figure()
    plt.plot(tableResultsSorted[:,4],'*',label='BQCD')
    plt.ylabel('BQCD')
    plt.legend()
    plt.show()
    
elif plotSelection in ['plotSorted10']:
#    glassA_idx  ,glassB_idx,  thetaA_final,  thetaB_final,  BDCQ_Final
    last = 3
    
    #Plot the prism angles
    plt.figure()    
    plt.plot(tableResultsSorted[0:last,2]/DEG2RAD,'*',label='theta A')
    plt.plot(tableResultsSorted[0:last,3]/DEG2RAD,'*',label='theta B')
    plt.ylabel('theta [°]')
    plt.legend()
    plt.show()
    
    #Plot the BDCQ value
    plt.figure()
    plt.plot(tableResultsSorted[0:last,4],'*',label='BQCD')
    plt.ylabel('BQCD')
    plt.legend()
    plt.show()
    
    for combIdx in range(last):
        print('Combination ', combIdx, ' : BDCQ = ',tableResultsSorted[combIdx,4],' | idx :', tableResultsSorted[combIdx,0], '/',tableResultsSorted[combIdx,1],
              ' | Glasses : ' , SellmeierCoefficients.GlassName[tableResultsSorted[combIdx,0]],
                '/' , SellmeierCoefficients.GlassName[tableResultsSorted[combIdx,1]])
        
    #Plot the refractive index wrt lambda for the best combination of glasse
    lambda_um = 1e6*np.linspace(lambda_wave[0],lambda_wave[-1],100)
    lambda_um = 1e6*np.linspace(300e-9,2400e-9,100)
    idxA = tableResultsSorted[0,0]
    idxB = tableResultsSorted[0,1]
    n_A = sellmeierRefractiveIndex(lambda_um,SellmeierCoefficients[int(idxA):int(idxA+1)])
    n_B = sellmeierRefractiveIndex(lambda_um,SellmeierCoefficients[int(idxB):int(idxB+1)])

    plt.figure()
    plt.plot(lambda_um, n_A,'+',label='n_A')
    plt.plot(lambda_um, n_B,'+',label='n_B')
    plt.ylabel('n(\lambda)')
    plt.legend()
    plt.show()


#tableResults[0:5,0]
#Out[10]: array([507., 507., 227., 105., 111.])
#
#tableResults[0:5,1]
#Out[11]: array([194., 177., 505., 237., 219.])
#
#SellmeierCoefficients.GlassName[507]
#Out[12]: 'S-NBH57'
#
#SellmeierCoefficients.GlassName[194]
#Out[13]: 'BAM21'





