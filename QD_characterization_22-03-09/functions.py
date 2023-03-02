# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 15:41:59 2022

@author: lucas_000
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

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

def maxima_reference(df_buffer_OD_1p0, number_maxima, cl_run, structure_name): 

    maxima = np.arange(0, number_maxima, 1)


    df = find_maxima(df_buffer_OD_1p0)


    df = sort_int(df)


    df_maxima_OD_1p0 = df[0:len(maxima)]
    df_maxima_OD_1p0.reset_index(level=0, inplace=True)

    
    maxima_frame = sort_wavelength(df_maxima_OD_1p0)
    maxima_frame.reset_index(level=0, inplace=True, drop=True)


    return maxima_frame, maxima    
#%%
def ps_maxima_matrix(background_frame, maxima_frame, power_series, directory, file_name, maxima):
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
        df_base[str(OD)] = df_base[str(OD)] - background_frame['Int']
        
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

    
    return ps_maxima_matrix

#%%
def pol_maxima_matrix(background_frame, maxima_frame, pol_series, directory, file_name, maxima):
    
    '''
    Gets the maxima for the identified peaks of 'maxima_frame' for varying polarizations
    '''
    pol_maxima_matrix = pd.DataFrame({'Wavelength': maxima_frame['Wavelength']})

    for pol in pol_series:
        pol_angle = (pol-1)*5
        real_file_name = file_name+str(pol)
        file_path =  os.path.join(directory, real_file_name+'.csv')
     
        df_base = pd.read_csv(file_path,
                       sep='\t', header = 0)
        
        df_base = df_base[['Wavelength', 'Intensity']]
        df_base.columns = ['Wavelength', str(pol_angle)]
        df_base[str(pol_angle)] = df_base[str(pol_angle)]-background_frame['Int']
        
        ps_pol = pd.DataFrame({'Wavelength': [0],
                              str(pol_angle): 0.0})
        for i in maxima:
            ix = int(maxima_frame.iloc[i]['index'])
            ps_pol_running = pd.DataFrame({'Wavelength': [df_base.iloc[ix]['Wavelength']],
                                          str(pol_angle): df_base.iloc[ix][str(pol_angle)]})
            ps_pol = pd.concat([ps_pol, ps_pol_running])
        ps_pol = ps_pol.reset_index(drop=1)
        ps_pol = ps_pol.iloc[1:]
        ps_pol = ps_pol.reset_index(drop=1)
    
        
        pol_maxima_matrix_running = pd.DataFrame({str(pol_angle): ps_pol[str(pol_angle)]})
        pol_maxima_matrix = pd.concat([pol_maxima_matrix, pol_maxima_matrix_running], axis=1)
    
    return pol_maxima_matrix
#%%
def pol_spectra_matrix(df_pol_0, pol_series, pol_step, directory, file_name, cl_run, structure_name):
    
    '''
    Gets the full spectra of the polarization series
    '''
    
    pol_spectra_matrix = pd.DataFrame({'Wavelength': df_pol_0['Wavelength']})

    for pol in pol_series:
        pol_angle = (pol-1)*pol_step
        real_file_name = file_name+str(pol)
        file_path =  os.path.join(directory, real_file_name+'.csv')
     
        df_base = pd.read_csv(file_path,
                       sep='\t', header = 0)
        
        df_base = df_base[['Wavelength', 'Intensity']]
        df_base.columns = ['Wavelength', str(pol_angle)]
        
        pol_spectra_matrix_running = pd.DataFrame({str(pol_angle): df_base[str(pol_angle)]})
        pol_spectra_matrix = pd.concat([pol_spectra_matrix, pol_spectra_matrix_running], axis=1)
        
    pol_zoom_plot = pol_spectra_matrix.drop(['Wavelength'], axis=1)
    wavelength = pol_spectra_matrix['Wavelength']
    angles = (pol_series-1)*pol_step
    
    fig, ax = plt.subplots(figsize=(10,4))
    e = [min(wavelength), max(wavelength), max(angles), min(angles)]
    cax = ax.matshow(pol_zoom_plot.transpose(), aspect="auto",norm=LogNorm(vmin=500, vmax=pol_zoom_plot.max().max()), extent=e, cmap='inferno')
    #cbar = fig.colorbar(cax, label="Intensity (arb. units)")

    ax.set_ylabel("Angle (°)")
    ax.set_xlabel("Wavelength (nm)")
    ax.set_xlim(min(wavelength), max(wavelength))
    ax.xaxis.set_ticks_position("bottom")
    plt.title('Pol Series '+str(cl_run)+' '+str(structure_name)+'.')
    
    plt.savefig(os.path.join(directory, str(cl_run)+'_'+str(structure_name))+'PolSeries_Full.png')
        
    return pol_spectra_matrix

#%%
    
def plot_pol_spectra_zoom(maxima, maxima_frame, pol_spectra_matrix, pol_series, pol_step, directory, cl_run, structure_name):
    '''
    Plots a heatmap around the identified peaks from the pol series spectral data
    '''
    for i in maxima:
        
        max_index = maxima_frame.iloc[i]['index']
        lower_index = int(max_index) - 10
        upper_index = int(max_index) + 10
        
        if lower_index < 0:
            lower_index = int(0)
        if upper_index > 1339:
            upper_index = int(1339)
        
        len_index = upper_index - lower_index
        
        pol_spectra_zoom = pol_spectra_matrix.iloc[lower_index:upper_index]
        
        pol_spectra_zoom.to_csv(os.path.join(directory, str(cl_run)+'_'+str(structure_name)+'PolSeries_Maximum '+str(i+1)+'.txt'),
                                index=None,
                                sep='\t',
                                mode='w')
        
        pol_zoom_plot = pol_spectra_zoom.drop(['Wavelength'], axis=1)
        wavelength = pol_spectra_zoom.iloc[len_index-len_index:len_index-1]['Wavelength']
        angles = (pol_series-1)*pol_step
        
        fig, ax = plt.subplots(figsize=(5,4))
        e = [min(wavelength), max(wavelength), max(angles), min(angles)]
        cax = ax.matshow(pol_zoom_plot.transpose(), aspect="auto",norm=LogNorm(vmin=500, vmax=pol_zoom_plot.max().max()), extent=e, cmap='inferno')
        #cbar = fig.colorbar(cax, label="Intensity (arb. units)")
    
        ax.set_ylabel("Angle (°)")
        ax.set_xlabel("Wavelength (nm)")
        ax.set_xlim(min(wavelength), max(wavelength))
        ax.xaxis.set_ticks_position("bottom")
        plt.title('Pol Series '+str(cl_run)+' '+str(structure_name)+' at Maximum '+str(i+1)+'.')
        
        plt.savefig(os.path.join(directory, str(cl_run)+'_'+str(structure_name))+'PolSeries_Maximum '+str(i+1)+'.png')
       
    return pol_spectra_zoom
    
#%%
def transpose_power_matrix(power_series, ps_maxima_matrix):
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
