# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 15:38:23 2022

@author: lucas_000
"""

'''
Extract Data from Pol Series Measurements
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from functions import maxima_reference, plot_line, plot_line_annotate, pol_maxima_matrix

#%%
root = r'C:\Users\lucas_000\tubCloud\Promotion\Proben_Messdaten'
sample_folder = 'np7654_7655_900nm_FlipChip_CBG'
pl_or_cl = 'PL'
mes_date = '2022-03-05'
mes_type = 'Polarization'
mes_type2 = 'pol'
cl_run = 'CL02'
structure_name = 'qd3' #map6p1

wafer = 'np7654'
opo_wl = '780nmOPO'
integ_time = '500ms'
grating = '1200g' 

pol_series = np.arange(1, 74, 1)

directory = os.path.join(root,
                         sample_folder,
                         pl_or_cl,
                         mes_date,
                         mes_type,
                         cl_run,
                         structure_name)
file_name = str(wafer)+'_'+str(cl_run)+'_'+str(structure_name)+'_'+str(opo_wl)+'-'+str(integ_time)+'-'+str(grating)+'_'+str(mes_type2)+'_  '

    
#%%
'''
Establish Reference Spectrum at OD 1.0 Power value
'''
   
pol_0 = 1
real_file_name = file_name+str(pol_0)
file_path =  os.path.join(directory, real_file_name+'.csv')
df_pol_0 = pd.read_csv(file_path,
                   sep='\t', header = 0)

    
df_pol_0 = df_pol_0[['Wavelength', 'Intensity']]
df_pol_0.columns = ["Wavelength", "Int"] 

num_maxima = 10     # number of expected maxima (ordered by intensity)

plot_line(df_pol_0['Wavelength'],
          df_pol_0['Int'],
          'Wavelength (nm)',
          'Intensity (a.u.)',
          'Pol 0° Raw Data '+str(cl_run)+' '+str(structure_name))

maxima_reference_res = maxima_reference(df_pol_0, num_maxima, cl_run, structure_name)
maxima_frame = maxima_reference_res[0]
maxima = maxima_reference_res[1]

plot_line_annotate(df_pol_0['Wavelength'],
                   df_pol_0['Int'],
                   'Wavelength (nm)',
                   'Intensity (a.u.)',
                   'Pol 0° Identified Maxima '+str(cl_run)+' '+str(structure_name),
                   maxima,
                   maxima_frame,
                   os.path.join(directory, str(cl_run)+'_'+str(structure_name)+'_Identified_Peaks.png'))

pol_maxima_matrix = pol_maxima_matrix(maxima_frame, 
                                     pol_series,
                                     directory,
                                     file_name,
                                     maxima)

#%%



pol_spectra_matrix = pd.DataFrame({'Wavelength': df_pol_0['Wavelength']})

for pol in pol_series:
    pol_angle = (pol-1)*5
    real_file_name = file_name+str(pol)
    file_path =  os.path.join(directory, real_file_name+'.csv')
 
    df_base = pd.read_csv(file_path,
                   sep='\t', header = 0)
    
    df_base = df_base[['Wavelength', 'Intensity']]
    df_base.columns = ['Wavelength', str(pol_angle)]
    ps_pol = pd.DataFrame({'Wavelength': [0],
                          str(pol_angle): 0.0})
    for i in maxima:
        ix = int(maxima_frame.iloc[i]['index'])
        ps_pol_running = pd.DataFrame({'Wavelength': [df_base.iloc['Wavelength']],
                                      str(pol_angle): df_base.iloc[str(pol_angle)]})
        ps_pol = pd.concat([ps_pol, ps_pol_running])
    ps_pol = ps_pol.reset_index(drop=1)
    ps_pol = ps_pol.iloc[1:]
    ps_pol = ps_pol.reset_index(drop=1)

    
    pol_maxima_matrix_running = pd.DataFrame({str(pol_angle): ps_pol[str(pol_angle)]})
    pol_maxima_matrix = pd.concat([pol_maxima_matrix, pol_maxima_matrix_running], axis=1)


