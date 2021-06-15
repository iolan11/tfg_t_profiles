# -*- coding: utf-8 -*-
"""
Created on Wed May  5 16:40:47 2021

@author: iolan
"""
import numpy as np
import math as ma

def Ziagos_func(K, z_exp, l_aq, bg_grad, palpha, pkappa):
      
    x = K[0]
    t = K[1]
    vf = K[2]
    kr = 2.65
    ta = 8.75

    alpha = palpha / vf
    kappa = (pkappa  * kr) / 4
    
    t2 = t - x/vf - (alpha * l_aq * x)/(3 * kappa)
    
    T = []
    
    for i in range(0, len(z_exp)):
        
        if z_exp[i] < l_aq:
            T.append(calc_T1(l_aq, z_exp[i], alpha, x, kappa, t2, bg_grad, ta))
        elif z_exp[i] > l_aq:
            T.append(calc_T2(l_aq, z_exp[i], alpha, x, kappa, t2, bg_grad, ta))
        elif z_exp[i] == l_aq:
            T.append(calc_T3(l_aq, alpha, x, kappa, t2, bg_grad, ta))

    return T

def calc_T1(l_aq, z_exp, alpha, x, kappa, t2, bg_grad, ta):
    
    sum_aux = calcSumT1(l_aq, z_exp, alpha, x, kappa, t2)
    T1_n = (ta * np.exp(-(alpha * x)/l_aq) * sum_aux)
    T1_grad = z_exp * bg_grad
    print(T1_n + T1_grad)
    return (T1_n + T1_grad)

def calc_T2(l_aq, z_exp, alpha, x, kappa, t2, bg_grad, ta):
    
    arg_aux = (alpha * x + z_exp - l_aq)/((4*kappa*t2)**(1/2))
    T2_n = (ta * np.exp(-(alpha * x)/l_aq) * ma.erfc(arg_aux))
    T2_grad = z_exp * bg_grad
    print(T2_n + T2_grad)
    return (T2_n + T2_grad)

def calc_T3(l_aq, alpha, x, kappa, t2, bg_grad, ta):
    
    arg_aux = ((alpha * x)/(4*kappa*t2)**(1/2))
    T3_n = (ta * np.exp(-(alpha * x)/l_aq) * ma.erfc(arg_aux))
    T3_grad = l_aq * bg_grad
    print(T3_n + T3_grad)
    return (T3_n + T3_grad)
    
    
def calcSumT1(l_aq, z, alpha, x, kappa, t2):
       
    suma = 0
    n = 1000
    
    for i in range(0, n):
        
        arg1 = ((2*i + 1)*l_aq - z + alpha*x)/((4*kappa*t2)**(1/2))      
        sum_add = ma.erfc(arg1)        
        
        arg2 = ((2*i + 1)*l_aq + z + alpha*x)/((4*kappa*t2)**(1/2))      
        sum_res = ma.erfc(arg2)
        
        suma = suma + sum_add - sum_res
                
    
    return suma
