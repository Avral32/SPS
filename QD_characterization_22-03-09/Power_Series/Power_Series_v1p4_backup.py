# -*- coding: utf-8 -*-
"""
Created on Fri Feb 25 14:40:17 2022

@author: lucas_000
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from functions import maxima_reference, plot_line_annotate, ps_maxima_matrix, transpose_matrix

#%%
root = r'C:\Users\lucas_000\tubCloud\Promotion\Proben_Messdaten'
sample_folder = 'np7654_7655_900nm_FlipChip_CBG'
pl_or_cl = 'PL'
mes_date = '2022-02-22'
mes_type = 'Power Series'
mes_type2 = 'power'
cl_run = 'CL02'
structure_name = 'map10p1' #map6p1

wafer = 'np7654'
opo_wl = '780nmOPO'
integ_time = '100ms'
grating = '1200g' 

power_series = np.arange(1, 31, 1)

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
   
power_OD_1p0 = 21
real_file_name = file_name+str(power_OD_1p0)
file_path =  os.path.join(directory, real_file_name+'.csv')
df_OD_1p0 = pd.read_csv(file_path,
                   sep='\t', header = 0)

    
df_OD_1p0 = df_OD_1p0[['Wavelength', 'Intensity']]
df_OD_1p0.columns = ["Wavelength", "Int"] 

num_maxima = 10     # number of expected maxima (ordered by intensity)

maxima_reference_res = maxima_reference(df_OD_1p0, num_maxima, cl_run, structure_name)
maxima_frame = maxima_reference_res[0]
maxima = maxima_reference_res[1]

plot_line_annotate(df_OD_1p0['Wavelength'],
                   df_OD_1p0['Int'],
                   'Wavelength (nm)',
                   'Intensity (a.u.)',
                   'OD 1.0 Identified Maxima '+str(cl_run)+' '+str(structure_name),
                   maxima,
                   maxima_frame,
                   os.path.join(directory, str(cl_run)+'_'+str(structure_name)+'_Identified_Peaks.png'))

ps_maxima_matrix = ps_maxima_matrix(maxima_frame, 
                                    power_series,
                                    directory,
                                    file_name,
                                    maxima)

WLs = ps_maxima_matrix['Wavelength']

ps_transpose_res = transpose_matrix(power_series,
                                ps_maxima_matrix)

ps_transpose = ps_transpose_res[0]
df_OD = ps_transpose_res[1]

print('')
ps_transpose.to_csv(os.path.join(directory, str(cl_run)+'_'+str(structure_name)+'_PowerSeries.txt'), index=None, sep='\t', mode='w')
print('Transpose Power Series Data saved! ')
print('')

''' 
Plot the extracted peak maxima intensities from the power series
''' 

plt.figure()
plt.xlabel('OD')
plt.ylabel('Peak Intensity (from raw data) a.u.')
i = 1
for WL in WLs: 
    
    ax = plt.gca()
    x_data = -df_OD
    y_data = ps_transpose[str(WL)]
    plt.scatter(x_data, y_data, lw=1, label=str(i))
    ax.set_yscale('log')
    plt.title('Power Series '+str(cl_run)+' '+str(structure_name))
    ax.legend(markerscale = 0.7, fontsize='small')
    plt.savefig(os.path.join(directory, str(cl_run)+'_'+str(structure_name)+'_PowerSeries.png'))
    plt.show
    i = i+1
