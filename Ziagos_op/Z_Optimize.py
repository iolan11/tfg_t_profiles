# -*- coding: utf-8 -*-
"""
Created on Wed May  5 16:40:34 2021

@author: iolan
"""
import scipy.optimize as optimize
import FileHandler as fh
import numpy as np
import matplotlib.pyplot as plt
import Z_Func as ZF

"""1. Recuperamos datos fichero de entrada"""

n_file = "S6_opt"
z_exp,T_exp = fh.readMyFile('data/' + n_file + '.csv')
for i in range(0,len(z_exp)):
        
    z_exp[i] = float(z_exp[i].replace(",", "."))    
    T_exp[i] = float(T_exp[i].replace(",", "."))
        
z_exp = np.array(z_exp, dtype=float)
T_exp = np.array(T_exp, dtype=float)

#Estandarizar los datos experimentales:
z_exp = z_exp - z_exp[0]
T_exp = T_exp - T_exp[0]

p1 = plt.figure()
simbol = u"\u00b0"
plt.xlabel('T('+simbol+'C)')
plt.ylabel('Depth(m)')
plt.gca().invert_yaxis()
plt.grid(color='grey')
plt.ylim(400, z_exp[0])
plt.xlim(T_exp[0], 22.5)
plt.plot(T_exp, z_exp, 'o', label="Experimental data", color="green")

#Ks = [x, t, vf, kr, ta]
#z_exp
l_aq = 90 #m
bg_grad = 0.0535 #ºC/m
#bg_grad = (z_exp[-1] - z_exp[0]) / (T_exp[-1] - T_exp[0])
print(bg_grad)
alpha = 1.26
kappa = 25.24

K0 = np.array([0, 1000, 1, 2, 20]) #valores iniciales
K_bounds = ((0, 1000), (1000, 50000), (1, 10), (2, 3.5), (20, 120))

K0 = np.array([50, 2000, 1, 2, 40]) #valores iniciales
K_bounds = ((50, 500), (2000, 10000), (1, 10), (2, 3.4), (40, 120))

K0 = np.array([50, 2000, 1, 2, 0]) #valores iniciales
K_bounds = ((50, 500), (2000, 10000), (1, 5), (2, 3.4), (0, 15))

K0 = np.array([50, 2000, 1, 2, 0]) #valores iniciales
K_bounds = ((50, 500), (2000, 5000), (1, 5), (2, 3.2), (0, 15))

K0 = np.array([50, 1000, 0.1]) #valores iniciales
K_bounds = ((50, 500), (1000, 5000), (0.1, 10))

K0 = np.array([50, 1000, 1]) #valores iniciales
K_bounds = ((50, 500), (1000, 500000), (1, 10))

K0 = np.array([50, 1000, 1]) #valores iniciales
K_bounds = ((50, 1000), (1000, 5000000), (1, 10))

K0 = np.array([50, 1000, 1]) #valores iniciales
K_bounds = ((50, 1000), (1000, 5e10), (1, 10))

K0 = np.array([50, 1000, 1]) #valores iniciales
K_bounds = ((50, 1000), (1000, 5e12), (1, 10))

kr = 2.65
ta = 8.75

def func_objetivo(Ks):  
   
    print(Ks)
    K = np.array([Ks[0], Ks[1], Ks[2], kr, ta])
    print(K)
    T = ZF.Ziagos_func(K, z_exp, l_aq, bg_grad, alpha, kappa)
    error = sum((T - T_exp)**2)
    print(Ks)
    return error

"""Optimizamos"""
result = optimize.minimize(func_objetivo, K0,  method='Powell', bounds=K_bounds, tol=1e-10)

print(result)
K = np.array([result.x[0], result.x[1], result.x[2], kr, ta])
T_result= ZF.Ziagos_func(result.x, z_exp, l_aq, bg_grad, alpha, kappa)
print(T_result)
plt.plot(T_result, z_exp, '-.', label="Adjusted model", color="red")
plt.legend(bbox_to_anchor=(1,1.02), ncol=1, fontsize='x-small')
plt.savefig('adj' + '.png')