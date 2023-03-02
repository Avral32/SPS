# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 15:51:22 2022

@author: lucas_000
"""

#%%
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

'''
Extract power series data for the whole series
for a given map
'''

cbg = '10p1'
cl_run = 'CL02' 
   
directory = r'C:\Users\lucas_000\tubCloud\Promotion\Bachelor_Master_Studenten\Rodrigo\Messungen\2022-02-22\Power Series'
file_name = 'np7654_'+str(cl_run)+'_map'+str(cbg)+'_780nmOPO-100ms-1200g_power_  '

power_OD_1p0 = 21
real_file_name = file_name+str(power_OD_1p0)
file_path =  os.path.join(directory, cbg, real_file_name+'.csv')
df_buffer_OD_1p0 = pd.read_csv(file_path,
                   sep='\t', header = 0)
df_buffer_OD_1p0 = df_buffer_OD_1p0[['Wavelength', 'Intensity']]
df_buffer_OD_1p0.columns = ["Wavelength", "Int"] 


power_series = np.arange(1,31,1)

ps_matrix = pd.DataFrame({'Wavelength': df_buffer_OD_1p0['Wavelength']})

for power in power_series:
    OD = 3.1-power*0.1
    OD = round(OD, 1)
    real_file_name = file_name+str(power)
    file_path =  os.path.join(directory, cbg, real_file_name+'.csv')
    print('')
    print('The directory is: ')
    print(directory)
    print('')
    print('The file_name is: ')
    print(real_file_name)
    df_base = pd.read_csv(file_path,
                   sep='\t', header = 0)

    df_base = df_base[['Wavelength', 'Intensity']]
    df_base.columns = ["Wavelength", str(OD)]
#    print(df_base)
    
    wl = df_base['Wavelength']
    Int = df_base[str(OD)]
    ps_matrix_running = df_base[str(OD)]
    print(ps_matrix_running)
    ps_matrix = pd.concat([ps_matrix, ps_matrix_running], axis=1)
    
    plt.figure()
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Intensity (a.u.)')
    plt.plot(wl, Int, lw=1)
    plt.title('Raw Data - OD: '+str(OD))
    plt.show
    
print(ps_matrix)
