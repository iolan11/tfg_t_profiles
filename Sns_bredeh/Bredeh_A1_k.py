# -*- coding: utf-8 -*-
"""
BREDEHOEFT AND PAPADOPULOS SOLUTION
"""
import matplotlib.pyplot as plt
import numpy as np

c_t0 = 0 #ºC

    
z_exp = np.linspace(0, 150, 20)
#Constant values
c_tL = 50   #ºC
c_cond_all = [2.2, 2.8, 3.3] #J/s m K (conductivitat tèrmica del sòlid)
c_cap = 4180000 #J/ K m^3 (calor específica del fluid)
c_L = 300 #m (longitud total)
vz = 0.3 #m/yr
c_gra = 0.03
    
T_gra = []
t_cnt = 14.5
    
for k in range(0,len(z_exp)):
    T_gra.append((z_exp[k] * c_gra) + t_cnt)
    
simbol = u"\u00b0"
plt.figure(dpi=1200)
plt.xlabel('T (' + simbol + 'C)', fontsize = 15)
plt.ylabel('Depth (m)', fontsize = 15)
plt.gca().invert_yaxis()
plt.grid(color='grey')
plt.xlim(0, 50)
plt.ylim(160, 0)
plt.yticks(np.arange(0, 180, step=40))
plt.text(45, 20, "b)",fontsize=14,fontweight='bold')

vz = -vz/31556926

for i in range(0, len(c_cond_all)):
    
   
    c_cond = c_cond_all[i]
    
    beta = (c_cap * vz * c_L) / c_cond
    
    fB_prt1 = np.exp((beta*z_exp)/c_L) - 1
    fB_prt2 = np.exp(beta) - 1
    fB = fB_prt1/fB_prt2
    T = fB * (c_tL - c_t0) + c_t0
    
    str_legend = str(c_cond) + " W/m"+simbol+"C"
    p = plt.plot(T, z_exp, '-', label = str_legend, linewidth = 2)
    
plt.legend(ncol=1, fontsize=10, loc = "lower left")
plt.savefig('k.png')