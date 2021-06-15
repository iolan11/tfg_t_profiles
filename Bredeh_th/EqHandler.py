# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 13:30:33 2021

@author: leoel
"""
import numpy as np


def calcTz(beta, depth, L, T0, TL):
    
    Tz = []
    
    for i in range(0,len(depth)):
               
        Tz.append(funcBeta(beta, depth[i], L) * (TL - T0) + T0)
    
    return Tz

 
def funcBeta(beta, z, L):
    
    fB_prt1 = np.exp((beta*z)/L) - 1
    fB_prt2 = np.exp(beta) - 1
    
    fB = fB_prt1/fB_prt2
    
    return fB

def calcBeta(heat_cap, vz, L, T_cond):
    
    beta = (heat_cap * vz * L) / T_cond
    
    return beta

def calcfBetaZL(zL, beta):
    
    fBeta = []
    
    for i in range(0, len(zL)):
        
        fBeta.append((np.exp(beta*(zL[i])) - 1)/(np.exp(beta) - 1))
          
    return fBeta

def vzByBeta(beta, c_cond, c_cap, L):
    
    vz = (beta * c_cond) / (c_cap * L)
          
    return vz

    
    