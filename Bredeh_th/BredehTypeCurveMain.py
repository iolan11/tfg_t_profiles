# -*- coding: utf-8 -*-
"""
BREDEHOEFT AND PAPADOPULOS SOLUTION: Z/L(f) comparando curvas
"""
import FileHandler as fh
import matplotlib.pyplot as plt
import numpy as np
import EqHandler as eq
from sklearn.metrics import mean_squared_error


"""Read file with experimental data - csv 2 columns"""
z_exp,T_exp = fh.readMyFile('data/S2_exp.csv')


for i in range(0,len(z_exp)):
    z_exp[i] = float(z_exp[i].replace(",", "."))
    T_exp[i] = float(T_exp[i].replace(",", "."))
  
z_exp = np.array(z_exp, dtype=float)
T_exp = np.array(T_exp, dtype=float)
  
#Constant values
c_L0 = 0 #m
c_LL = z_exp[-1] #m (depending on the well)
c_L = c_LL - c_L0 #m (longitud total)
c_t0 = 14 #ºC
c_tL = T_exp[-1]   #ºC
c_cond = 2.5 #J/s m K (conductivitat tèrmica del sòlid)
c_cap = 4180000 #J/ K m^3 (calor específica del fluid)

zL_exp = z_exp/abs(c_L)
x_exp  = (T_exp - c_t0) / (c_tL - c_t0)

p = plt.plot(x_exp, zL_exp, 'o')
plt.xlabel('f')
plt.ylabel('z/L')
plt.gca().invert_yaxis()
    
"""Por cada valor posible de z/L calculamos la curva f(beta, -z/L) y 
miramos los errores de las curvas y cogemos el mínimo"""

zL_teo = zL_exp

#Calculamos beta_teo en dos partes porque nos tenemos que saltar el 0
beta_teo_1 = np.linspace(-30, -0.5, len(x_exp)//2)
beta_teo_2 = np.linspace(0.5, 30, len(x_exp)//2)
beta_teo = np.concatenate((beta_teo_1, beta_teo_2), axis=0)

d = {} #Diccionario con key: beta, valor: error

for i in range(0, len(beta_teo)):
    
    fBeta = eq.calcfBetaZL(zL_exp, beta_teo[i])
    d[beta_teo[i]] = np.sqrt(mean_squared_error(x_exp, fBeta))

#Buscamos el error mínimo
d_values = list(d.values())
d_keys = list(d.keys())
min_value = min(d_values)

beta_op = d_keys[d_values.index(min_value)] #beta más optima, con mínimo error
fBeta = eq.calcfBetaZL(zL_exp, beta_op)
print(beta_op)

vz = (beta_op * c_cond) / (c_cap * c_L)
print(vz)

p = plt.plot(fBeta, zL_exp, 'o')