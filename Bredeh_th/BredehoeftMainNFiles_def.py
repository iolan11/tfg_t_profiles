"""
BREDEHOEFT AND PAPADOPULOS SOLUTION - N Files
"""
import FileHandler as fh
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import mean_squared_error, r2_score
import pandas as pd

"""Read file with datafiles names"""
n_files, n_t0 = fh.readMyFile('data/Files_List.csv')

print(n_files)

for k in range(0,len(n_t0)):
        
        n_t0[k] = float(n_t0[k].replace(",", "."))    
        
n_t0 = np.array(n_t0, dtype=float)

tabla = pd.DataFrame(columns=["Well", "vz (m/yr)", "error"])

t_cnt = 14.5

for j in range(0, len(n_files)):
        
    """Read file with experimental data - csv 2 columns"""
    z_exp,T_exp = fh.readMyFile('data/' + n_files[j] + '.csv')
    
    c_t0 = n_t0[j]
    
    for i in range(0,len(z_exp)):
        
        z_exp[i] = float(z_exp[i].replace(",", "."))    
        T_exp[i] = float(T_exp[i].replace(",", "."))
        
    z_exp = np.array(z_exp, dtype=float)
    T_exp = np.array(T_exp, dtype=float)
    
    """Show experimental data"""
    plt.figure(j)
    p = plt.plot(T_exp, z_exp, 'o', label="Experimental data", color="green")
    plt.xlabel('T(ºC)')
    plt.ylabel('Depth(m)')
    plt.gca().invert_yaxis()
    
    plt.grid(color='grey')
    
    """Make gradient"""
    c_gra = 0.03
    
    T_gra = []
    
    for k in range(0,len(z_exp)):
        T_gra.append((z_exp[k] * c_gra) + t_cnt)
        #T_gra.append((z_exp[k] * c_gra))
    
    p = plt.plot(T_gra, z_exp, '-', label="Background gradient", color="gray")
    
    """Make theorical data"""
    z_teo = z_exp 
    
    #Constant values
    c_tL = T_exp[-1]   #ºC
    c_L0 = 0 #m
    c_LL = z_exp[-1] #m (depending on the well)
    c_cond = 2.65 #J/s m K (conductivitat tèrmica del sòlid)
    c_cap = 4180000 #J/ K m^3 (calor específica del fluid)
    c_L = c_LL - c_L0 #m (longitud total)
    
    #Iremos dando diferentes valores a la velocidad. Por cada valor, calcularemos
    #T_teo y en función del error entre datos experimentales y teoricos
    #seguiremos dando valores a la velocidad o nos quedaremos con ese
    
    vz_pre = 0
    n_vz = True
    
    while n_vz:
        
        vz_pre = vz_pre - 0.00000000001
         
        # beta_pre = eq.calcBeta(c_cap, vz_pre, abs(c_L), c_cond)
        beta_pre = (c_cap * vz_pre * c_L) / c_cond
        
        # T_teo_pre = eq.calcTz(beta_pre, z_teo, c_L, c_t0, c_tL)
        fB_prt1 = np.exp((beta_pre*z_teo)/c_L) - 1
        fB_prt2 = np.exp(beta_pre) - 1
        fB = fB_prt1/fB_prt2
        T_teo_pre = fB * (c_tL - c_t0) + c_t0
    
        error_pre = np.sqrt(mean_squared_error(T_exp, T_teo_pre))
        r2_pre = r2_score(T_exp, T_teo_pre)
        
        
        vz_post = vz_pre - 0.00000000001
         
        # beta_post = eq.calcBeta(c_cap, vz_post, abs(c_L), c_cond)
        beta_post = (c_cap * vz_post * c_L) / c_cond
        
        # T_teo_post = eq.calcTz(beta_post, z_teo, c_L, c_t0, c_tL)   
        fB_prt1 = np.exp((beta_post*z_teo)/c_L) - 1
        fB_prt2 = np.exp(beta_post) - 1
        fB = fB_prt1/fB_prt2
        T_teo_post = fB * (c_tL - c_t0) + c_t0
        
        error_post = np.sqrt(mean_squared_error(T_exp, T_teo_post)) 
        r2_post = r2_score(T_exp, T_teo_post) 
        
        if error_post > error_pre: 
            
            n_vz = False
            
        
    print("vz " + n_files[j])
    print("%e" %vz_pre)
      
    fila = pd.DataFrame([[n_files[j], vz_pre*31556926, r2_pre]], columns=["Well", "vz (m/yr)", "error"])
    
    tabla = tabla.append(fila)
    
    fin = len(n_files[j])-4
    #if(r2_pre > 0.98):
    str_title = n_files[j][:fin] + " - R^2 = " + str(r2_pre)[0:7]
    #else:
       # str_title = n_files[j]
        
    plt.title(str_title)
    p = plt.plot(T_teo_pre, z_teo, '-.', label="Theoretical data", color="red")
    fig = plt.gcf()
    fig.legend()
    fig.savefig("out/" + n_files[j] + ".jpg")

csv_name = tabla.to_csv("data/well_vz.csv", index=False)

print(tabla)