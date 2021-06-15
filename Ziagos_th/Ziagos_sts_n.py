# -*- coding: utf-8 -*-
"""
ZIAGOS SOLUTION
"""
import matplotlib.pyplot as plt
import numpy as np
import math as ma

def calcSumT1(l_aq, z, alpha, x, kappa, t1):
    
    suma = 0
    n = 100000
    
    for i in range(0, n):
        
        arg1 = ((2*i + 1)*l_aq - z + 2*alpha*x)/((4*kappa*t1)**(1/2))      
        sum_add = ma.erfc(arg1)
        
        arg2 = ((2*i + 1)*l_aq + z + 2*alpha*x)/((4*kappa*t1)**(1/2))      
        sum_res = ma.erfc(arg2)
        
        suma = suma + sum_add - sum_res
        
        i = i + 1
    
    return suma

"""Long time aproximation"""
x_all = [0, 100, 250, 500, 1000]
t = 1200 #yr
l_aq = 100 #m
vf = 1 #m/yr
ta = 50 #ºC
kr = 4 #mcl/cm s ºC
bg_grad = 100 #ºC/km
a_th = 10 #m aquifer thickness
alpha = 1.26
ts = 0 #ºC surface temperature
kappa = 25.24 #todo difusivitat roca
z_max = 300

plt.xlabel('T(ºC)')
plt.ylabel('Depth(m)')
plt.gca().invert_yaxis()
plt.grid(color='grey')
plt.title("Short-time approximation")

for l in range(0, len(x_all)):
    
    x = x_all[l]
    
    t1 = t - x/vf
    
    """Càlcul de T1 va de 0 fins a l'aquifer"""
    z_T1 = np.linspace(0, l_aq, 20)
    T1 = []
    
    for i in range(0, len(z_T1)):
        
        sum_aux = calcSumT1(l_aq, z_T1[i], alpha, x, kappa, t1)
        T1.append(ta * sum_aux)
           
        
    """Càlcul de T2 va de l+a_th fins al final"""
    z_T2 = np.linspace(l_aq + a_th, z_max, 20)
    T2 = []
    
    for j in range(0, len(z_T2)):
        
        arg_aux = (2 * alpha * x + z_T2[j] - l_aq)/((4*kappa*t1)**(1/2))
        T2.append(ta * ma.erfc(arg_aux))
        
    
    """Càlcul de T3 a l """
    z_T3 = [100]
    T3 = []
            
    arg_aux = ((2 * alpha * x)/(4*kappa*t1)**(1/2))
    T3.append(ta * ma.erfc(arg_aux))
        
    T = T1 + T3 + T2
    z = z_T1.tolist() + z_T3 + z_T2.tolist()
    
    str_legend = str(x) + " m"
        
    p = plt.plot(T, z, '-', label = str_legend)
    
plt.legend()
