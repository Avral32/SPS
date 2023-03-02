# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 15:41:59 2022

@author: lucas_000
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#%%

def find_maxima(df):
    ''' Filters to the maximum intensity by comparing the 2nd nearest neighbours 
        Reduces the data to the found maxima
    '''
    df['max'] = df.Int[
    (df.Int.shift(1) < df.Int) & (df.Int.shift(-1) < df.Int) 
    & 
    (df.Int.shift(2) < df.Int) & (df.Int.shift(-2) < df.Int)
    ]
    
    # Reduce the wavelength vs intensity data to the found maxima
    df = df.dropna(how='any')
    df = df.loc[:, ['Wavelength', 'Int']]
    return df

#%%
    
def sort_int(df):
    ''' Sorts the Intensity data highest to lowest
    '''
     
    # Sorting the 'Int' values highest to lowest
    df = df.sort_values(by=['Int'], ascending=False)    
    
    return df

#%%
    
def re_index(df):
    ''' Re-indexes a WL-Int dataframe used for Maxima finding, 
    dropping the old index and reducing it to the WL and Int columns 
    '''
    df = df.reset_index(drop=0)
    df = df.loc[:, ['Wavelength', 'Int']]
    
    return df

#%%

def re_index_2(df):
    ''' Re-index any data frame and dropping the first row (typically a dummy row)
    '''
    df.reset_index(drop=0)
    df = df.iloc[1:]
    df.reset_index(drop=0)
    return df

#%%
    
def sort_wavelength(df):
    ''' Sorts a wavelength data frame shortest to longest 
    '''
    # Sorting the 'Wavelength' values longest to shortest
    df = df.sort_values(by=['Wavelength'], ascending=True)

    return df

#%%
    
def plot_line(x_data, y_data, xlabel, ylabel, title):
    ''' Plots a simple x-y line plot with lw 1
    passed values are x- and y-datapoints, x- and y-label (strings) 
    and title (string) of the plot
    '''
    
    plt.figure()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.plot(x_data, y_data, lw=1)
    plt.title(title)
    plt.show
    
#%%
    
def plot_line_annotate(x_data, y_data, xlabel, ylabel, title, num_annotates, df_annotate, file_path_name):
    ''' Simple x-y line plot with lw=1 including labels and title,
    and annotations at given x-y points of Dataframe 'annotates' 
    '''
    plt.figure()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.plot(x_data, y_data, lw=1)
    plt.title(title)
    for i in num_annotates:
        x_wl = df_annotate.iloc[i][1]
        y_int = df_annotate.iloc[i][2]
        x_offset = -0.2
        y_offset = df_annotate.max()[2]*0.01
        plt.annotate(str(i+1),
                     (x_wl, y_int),
                     (x_wl+x_offset, y_int+y_offset),
                     fontsize = 'xx-small')
    plt.savefig(file_path_name)
    plt.show

#%%

def maxima_reference(df_buffer_OD_1p0, number_maxima): 

    plot_line(df_buffer_OD_1p0['Wavelength'],
              df_buffer_OD_1p0['Int'],
              'Wavelength (nm)',
              'Intensity (a.u.)',
              'OD 1.0 Raw Data '+str(cl_run)+' '+str(structure_name))

    maxima = np.arange(0, number_maxima, 1)


    df = find_maxima(df_buffer_OD_1p0)


    df = sort_int(df)


    df_maxima_OD_1p0 = df[0:len(maxima)]
    df_maxima_OD_1p0.reset_index(level=0, inplace=True)

    
    maxima_frame = sort_wavelength(df_maxima_OD_1p0)
    maxima_frame.reset_index(level=0, inplace=True, drop=True)


    return maxima_frame, maxima    
#%%
def ps_maxima_matrix(maxima_frame, power_series, directory, file_name, maxima):
    '''
    Extract Maximum Intensity for the found maxima positions 
    and organize it into a matrix
    '''
    ps_maxima_matrix = pd.DataFrame({'Wavelength': maxima_frame['Wavelength']})

    
    for power in power_series:
        OD = 3.1-power*0.1
        OD = round(OD, 1)
        real_file_name = file_name+str(power)
        file_path =  os.path.join(directory, real_file_name+'.csv')
 
        df_base = pd.read_csv(file_path,
                       sep='\t', header = 0)
        
        df_base = df_base[['Wavelength', 'Intensity']]
        df_base.columns = ['Wavelength', str(OD)]
        ps_OD = pd.DataFrame({'Wavelength': [0],
                              str(OD): 0.0})
        for i in maxima:
            ix = int(maxima_frame.iloc[i]['index'])
            ps_OD_running = pd.DataFrame({'Wavelength': [df_base.iloc[ix]['Wavelength']],
                                          str(OD): df_base.iloc[ix][str(OD)]})
            ps_OD = pd.concat([ps_OD, ps_OD_running])
        ps_OD = ps_OD.reset_index(drop=1)
        ps_OD = ps_OD.iloc[1:]
        ps_OD = ps_OD.reset_index(drop=1)

        
        ps_maxima_matrix_running = pd.DataFrame({str(OD): ps_OD[str(OD)]})
        ps_maxima_matrix = pd.concat([ps_maxima_matrix, ps_maxima_matrix_running], axis=1)

    
    return(ps_maxima_matrix)

#%%
    
#%%
def transpose_matrix(power_series, ps_maxima_matrix):
    ''' 
    Transpose maxima_matrix for easier plotting
    '''
    ODs = 3.1-power_series*0.1
    
    WLs = ps_maxima_matrix['Wavelength']
    ps_transpose = pd.DataFrame({'OD': ODs})
    df_OD = ps_transpose.round(1)
    counter = np.arange(0, len(WLs), 1)
    counter_od = np.arange(0, len(ODs+1), 1)
    
    for i in counter:
        ps_OD_transpose = pd.DataFrame({'OD': [0.0],
                                        str(WLs[i]): [0]})
        for j in counter_od:
            ps_OD_transpose_running = pd.DataFrame({'OD': df_OD.iloc[j]['OD'],
                                                    str(WLs[i]): [ps_maxima_matrix.iloc[i][str(df_OD.iloc[j]['OD'])]]})
        
            ps_OD_transpose = pd.concat([ps_OD_transpose, ps_OD_transpose_running])
            
        ps_OD_transpose = ps_OD_transpose.reset_index(drop=1)
        ps_OD_transpose = ps_OD_transpose.iloc[1:]
        ps_OD_transpose = ps_OD_transpose.reset_index(drop=1)
        
        
        ps_transpose_running = pd.DataFrame({str(WLs[i]): ps_OD_transpose[str(WLs[i])]})
        
        ps_transpose = pd.concat([ps_transpose, ps_transpose_running], axis=1)

    
    return ps_transpose, df_OD

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

maxima_reference_res = maxima_reference(df_OD_1p0, num_maxima)
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
