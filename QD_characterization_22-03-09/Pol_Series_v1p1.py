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
from matplotlib.colors import LogNorm
import os
from functions import maxima_reference, plot_line, plot_line_annotate, pol_maxima_matrix, pol_spectra_matrix, plot_pol_spectra_zoom

#%%
root = r'C:\Users\lucas_000\tubCloud\Promotion\Proben_Messdaten'
sample_folder = 'np7654_7655_900nm_FlipChip_CBG'
pl_or_cl = 'PL'
mes_date = '2022-03-03'
mes_type = 'Polarization'
mes_type2 = 'pol'
cl_run = 'CL07'
structure_name = 'map7p2' #map6p1

wafer = 'np7654'
opo_wl = '780nmOPO'
integ_time = '100ms'
grating = '1200g' 

pol_series = np.arange(1, 74, 1)    # number of polarization files/measurements

pol_step = 5                        # Degree value of polarization step

directory = os.path.join(root,
                         sample_folder,
                         pl_or_cl,
                         mes_date,
                         mes_type,
                         cl_run,
                         structure_name)
file_name = str(wafer)+'_'+str(cl_run)+'_'+str(structure_name)+'_'+str(opo_wl)+'-'+str(integ_time)+'-'+str(grating)+'_'+str(mes_type2)+'_  '

mes_type_background = 'power_Test_OD1p0_Pol0deg_  1'

background_file_path = os.path.join(root,
                                    sample_folder,
                                    pl_or_cl)

background_file_name = 'Background_'+str(opo_wl)+'-'+str(integ_time)+'-'+str(grating)+'_'+str(mes_type_background)
    
#%%

'''
Get the Background Reference
'''


df_background = pd.read_csv(os.path.join(background_file_path,
                            background_file_name)+'.csv',
                            sep='\t',
                            header=0)
 
df_background = df_background[['Wavelength', 'Intensity']]
df_background.columns = ["Wavelength", "Int"]

#%%

'''
Establish Reference Spectrum at OD 1.0 Power value and Pol = 0°
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

pol_maxima_matrix = pol_maxima_matrix(df_background,
                                      maxima_frame, 
                                      pol_series,
                                      directory,
                                      file_name,
                                      maxima)

pol_spectra_matrix = pol_spectra_matrix(df_pol_0, 
                                        pol_series,
                                        pol_step,
                                        directory,
                                        file_name,
                                        cl_run,
                                        structure_name)
    
print('')
pol_spectra_matrix.to_csv(os.path.join(directory, str(cl_run)+'_'+str(structure_name)+'_PolSeries.txt'), index=None, sep='\t', mode='w')
print('Full Spectra of Polarization Series Data saved! ')
print('')

plot_pol_spectra_zoom(maxima,
                      maxima_frame,
                      pol_spectra_matrix,
                      pol_series,
                      pol_step,
                      directory,
                      cl_run,
                      structure_name)
    
print('Zoom-in Spectra of Polarization Series Data saved at respective Maxima! ')
print('')   
