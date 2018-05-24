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


alpha1, alpha0, V, h, pi, po,FP2PUPIMG, RPUPIMG,RTT= symbols('alpha1, alpha0, V, h, pi, po,FP2PUPIMG, RPUPIMG,RTT')
init_printing(use_unicode=True)

system = [alpha1 -alpha0+V*h, (FP2PUPIMG-po)*tan(alpha0)-h, 1/pi-1/po-V, h-RTT-pi*tan(alpha1),pi/po-RTT/RPUPIMG]
vars =  [alpha1, V, h, pi, po]
nonlinsolve(system,vars)