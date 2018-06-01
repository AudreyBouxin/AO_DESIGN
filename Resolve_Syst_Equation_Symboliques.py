# -*- coding: utf-8 -*-
"""
Created on Thu May 24 09:00:19 2018

@author: audrey.bouxin
"""

#FP2PUPIMG = 10289.737904
#RPUPIMG = 371.0345
#RTT = 5

from sympy import symbols
from sympy import tan
from sympy import *


#alpha1, alpha0, V, h, pi, po,FP2PUPIMG, RPUPIMG,RTT= symbols('alpha1, alpha0, V, h, pi, po,FP2PUPIMG, RPUPIMG,RTT')
#init_printing(use_unicode=True)
#
#system = [alpha1 -alpha0+V*h, (FP2PUPIMG-po)*tan(alpha0)-h, 1/pi-1/po-V, h-RTT-pi*tan(alpha1),pi/po-RTT/RPUPIMG]
#vars =  [alpha1, V, h, pi, po]
#nonlinsolve(system,vars)

alpha,DTT,DCCD,DLR,Pofoc,Pifoc,Pott,Pitt,fLR,dtt2foc= symbols('alpha,DTT,DCCD,DLR,Pofoc,Pifoc,Pott,Pitt,fLR,dtt2foc')
init_printing(use_unicode=True)

system = [DLR/DTT-Pofoc/dtt2foc, DCCD/DLR-(Pifoc-Pitt)/Pifoc,1/Pifoc-1/Pofoc-1/fLR,1/Pitt-1/Pott-1/fLR, Pott-Pofoc-dtt2foc,dtt2foc-DTT/tan(alpha)]
vars =  [DLR,Pofoc,Pifoc,Pott,Pitt,dtt2foc]
sol = nonlinsolve(system,vars)
sol = simplify(sol)
#print(sol)

alphaNUM = 1/60.26
DTTNUM=10
DCCDNUM = -0.792
fLRNUM = 60
DLR   = (sol.args[0][0]).evalf(subs={alpha:alphaNUM, DTT:DTTNUM, DCCD:DCCDNUM,fLR:fLRNUM})
Pofoc = (sol.args[0][1]).evalf(subs={alpha:alphaNUM, DTT:DTTNUM, DCCD:DCCDNUM,fLR:fLRNUM})
Pifoc = (sol.args[0][2]).evalf(subs={alpha:alphaNUM, DTT:DTTNUM, DCCD:DCCDNUM,fLR:fLRNUM})
Pott  = (sol.args[0][3]).evalf(subs={alpha:alphaNUM, DTT:DTTNUM, DCCD:DCCDNUM,fLR:fLRNUM})
Pitt  = (sol.args[0][4]).evalf(subs={alpha:alphaNUM, DTT:DTTNUM, DCCD:DCCDNUM,fLR:fLRNUM})
dtt2foc = (sol.args[0][5]).evalf(subs={alpha:alphaNUM, DTT:DTTNUM, DCCD:DCCDNUM,fLR:fLRNUM})

#print('DLR : ',DLR)
print('Pofoc : ',Pofoc)
print('Pifoc : ',Pifoc)
print('Pott : ',Pott)
print('Pitt : ',Pitt)
print('dtt2foc : ',dtt2foc)

#alpha=1/60
#DTT=10
#DCCD=-0.792
#fLR=60
#DLR = -DTT - fLR*tan(alpha) - DTT*fLR*tan(alpha)/DCCD
#Pofoc =  -DTT/tan(alpha) - fLR - DTT*fLR/DCCD
#Pifoc = fLR*(DCCD*DTT + DCCD*fLR*tan(alpha) + DTT*fLR*tan(alpha))/(DTT*(DCCD + fLR*tan(alpha)))
#Pitt = -fLR*(DCCD + DTT)/DCCD, fLR*(DCCD + DTT)/DTT
#dtt2foc = DTT/tan(alpha)
#print(Pofoc)




