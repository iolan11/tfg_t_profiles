# -*- coding: utf-8 -*-
"""
ZIAGOS SOLUTION
"""
import matplotlib.pyplot as plt
import numpy as np
import math as ma

def calcSumT1(l_aq, z, alpha, x, kappa, t2):
    
    suma = 0
    n = 100000
    
    for i in range(0, n):
        
        arg1 = ((2*i + 1)*l_aq - z + alpha*x)/((4*kappa*t2)**(1/2))      
        sum_add = ma.erfc(arg1)        
        
        arg2 = ((2*i + 1)*l_aq + z + alpha*x)/((4*kappa*t2)**(1/2))      
        sum_res = ma.erfc(arg2)
        
        suma = suma + sum_add - sum_res
                
    
    return suma

"""Long time aproximation"""
x_all = [0, 100, 250, 500, 1000]
t = 1500 #yr
l_aq = 100 #m
vf = 1 #m/yr
ta = 50 #ºC
kr = 4 #mcl/cm s ºC
bg_grad = 0.1 #ºC/m
a_th = 10 #m aquifer thickness
alpha = 1.26
ts = 0 #ºC surface temperature
kappa = 25.24 #todo difusivitat roca
z_max = 300

plt.xlabel('T(ºC)')
plt.ylabel('Depth(m)')
plt.gca().invert_yaxis()
plt.grid(color='grey')
plt.title("Long-time approximation")
plt.xlim(0, 100)

for l in range(0, len(x_all)):
    
    x = x_all[l]
    
    t2 = t - x/vf - (alpha * l_aq * x)/(3 * kappa)
    
    """Càlcul de T1 va de 0 fins a l'aquifer"""
    z_T1 = np.linspace(0, l_aq, 20)
    T1 = []
    
    for i in range(0, len(z_T1)):
        
        sum_aux = calcSumT1(l_aq, z_T1[i], alpha, x, kappa, t2)
        T1_n = (ta * np.exp(-(alpha * x)/l_aq) * sum_aux)
        T1_grad = z_T1[i] * bg_grad
        T1.append(T1_n + T1_grad)
        
       
    """Càlcul de T2 va de l+a_th fins al final"""
    z_T2 = np.linspace(l_aq, z_max, 20)
    T2 = []
    
    for j in range(0, len(z_T2)):
        
        arg_aux = (alpha * x + z_T2[j] - l_aq)/((4*kappa*t2)**(1/2))
        T2_n = (ta * np.exp(-(alpha * x)/l_aq) * ma.erfc(arg_aux))
        T2_grad = z_T2[j] * bg_grad
        T2.append(T2_n + T2_grad)
        
    
    """Càlcul de T3 a l"""
    z_T3 = [100]
    T3 = []
    
    arg_aux = ((alpha * x)/(4*kappa*t2)**(1/2))
    T3_n = (ta * np.exp(-(alpha * x)/l_aq) * ma.erfc(arg_aux))
    T3_grad = 100 * bg_grad
    T3.append(T3_n + T3_grad)
    
    T = T1 + T3 + T2
    z = z_T1.tolist() + z_T3 + z_T2.tolist()
    
    str_legend = str(x) + " m"
    
    p = plt.plot(T, z, '-', label = str_legend)
    
plt.legend()