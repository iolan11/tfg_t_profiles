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

"""Datos"""
ta_all = np.linspace(40, 120, 5) #ºC

kr = 4 #mcl/cm s ºC
vf = 1 #m/yr
x = 100 #m
t = 3000 #yr

l_aq = 100 #m
bg_grad = 0.1 #ºC/m
alpha = 1.26
kappa = 25.24 #todo difusivitat roca
z_max = 300

p1 = plt.figure()
plt.xlabel('T(ºC)')
plt.ylabel('Depth(m)')
plt.gca().invert_yaxis()
plt.grid(color='grey')
plt.tight_layout()

for l in range(0, len(ta_all)):
    
    ta = ta_all[l]
    
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
    z_T3 = [l_aq]
    T3 = []
    
    arg_aux = ((alpha * x)/(4*kappa*t2)**(1/2))
    T3_n = (ta * np.exp(-(alpha * x)/l_aq) * ma.erfc(arg_aux))
    T3_grad = l_aq * bg_grad
    T3.append(T3_n + T3_grad)
    
    T = T1 + T3 + T2
    z = z_T1.tolist() + z_T3 + z_T2.tolist()
    
    str_legend = str(ta) + " ºC"
    
    p = plt.plot(T, z, '-', label = str_legend)
    
plt.legend(loc='upper left',bbox_to_anchor=(1,1.02), ncol=1, fontsize='x-small')