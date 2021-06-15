# -*- coding: utf-8 -*-
"""
ZIAGOS SOLUTION
"""
import matplotlib.pyplot as plt
import numpy as np
import math as ma

def calcSumT1(l, z, alpha, x, kappa, t1):
    
    suma = 0
    n = 100000
    i = 0
    
    while i < n:
        
        arg1 = ((2*i + 1)*l - z + 2*alpha*x)/((4*kappa*t1)**(1/2))      
        sum_add = ma.erfc(arg1)
        
        arg2 = ((2*i + 1)*l + z + 2*alpha*x)/((4*kappa*t1)**(1/2))      
        sum_res = ma.erfc(arg2)
        
        suma = suma + sum_add - sum_res
        
        i = i + 1
    
    return suma

"""Long time aproximation"""
x = 0 #m
t = 1000 #yr
l = 100 #m
vf = 1 #m/yr
ta = 50 #ºC
kr = 4 #mcl/cm s ºC
bg_grad = 100 #ºC/km
a_th = 10 #m aquifer thickness
alpha = 1.26
ts = 0 #ºC surface temperature
kappa = 31.65 #todo difusivitat roca
z_max = 300

t1 = t - x/vf

"""Càlcul de T1 va de 0 fins a l'aquifer"""
z_T1 = np.linspace(0, l, 20)
T1 = []

for i in range(0, len(z_T1)):
    
    sum_aux = calcSumT1(l, z_T1[i], alpha, x, kappa, t1)
    T1.append(ta * sum_aux)
    
plt.xlabel('T(ºC)')
plt.ylabel('Depth(m)')
plt.gca().invert_yaxis()
plt.grid(color='grey')
plt.title("Short-time approximation")
    
"""Càlcul de T2 va de l+a_th fins al final"""
z_T2 = np.linspace(l + a_th, z_max, 20)
T2 = []

for j in range(0, len(z_T2)):
    
    arg_aux = (2 * alpha * x + z_T2[j] - l)/((4*kappa*t1)**(1/2))
    T2.append(ta * ma.erfc(arg_aux))
    

"""Càlcul de T3 va de l fins a l+a_th """
z_T3 = np.linspace(l, l + a_th, 20)
T3 = []

for k in range(0, len(z_T3)):
    
    arg_aux = ((2 * alpha * x)/(4*kappa*t1)**(1/2))
    T3.append(ta * ma.erfc(arg_aux))
    
T = T1 + T2 + T3
z = z_T1.tolist() + z_T2.tolist() + z_T3.tolist()
    
p = plt.plot(T, z, 'o-')
