# -*- coding: utf-8 -*-
"""
Created on Fri Dec 22 15:22:32 2023

@author: Gino Del Ferraro
"""

import numpy as np
import os
import pickle
# from utilities_ketamine_analysis_v8 import *
from utils_signal_processing import *
from utils_plotting import *
from utils_general import *
from scipy.signal import morlet, welch, convolve
import matplotlib.gridspec as gridspec
import time 

""" Upsample behaviour, i.e. speed, x, y, at 250 Hz from 100 Hz """

def upsample_behaviour(spk, speed, x, y, sess, spk_rate = 250, behav_rate = 100):
    
    # time length for speed and spike variables
    behav_T = np.linspace(0, len(speed[sess])/ behav_rate, len(speed[sess]))
    spk_T = np.linspace(0, spk[0].shape[1]/ spk_rate, spk[0].shape[1]) # all cells have same time length

    # interpolate speed variable based on spike 
    interpolator = interp1d(behav_T, speed[sess], kind = 'linear', fill_value="extrapolate")
    speed_upsampled = interpolator(spk_T)
    
    # interpolate x variable based on spike 
    interpolator = interp1d(behav_T, x[sess], kind = 'linear', fill_value="extrapolate")
    x_upsampled = interpolator(spk_T)
    
    # interpolate y variable based on spike 
    interpolator = interp1d(behav_T, y[sess], kind = 'linear', fill_value="extrapolate")
    y_upsampled = interpolator(spk_T)
    
    # check if upsample done correctly
    if (speed_upsampled.size - spk[0].shape[1]) !=0:
        sys.exit("Speed upsampled size and spike size are not the same!")
    print('speed upsampled shape {}, spike shape {}'.format(speed_upsampled.shape, spk[0].shape))
    
    return speed_upsampled, x_upsampled, y_upsampled



""" 
Compute instantaneous phase by using Morlet wavelet with n cycles, for a range of frequencies,
at a given sampling frequency for the signal for each minute and each channel of data 

Inputs: 
segment: LFP segment data for low speed
frequencies: range of frequencies for phase extraction (1D array)
num_cycle: n parameter in the Gaussian kernel for the number of cycles
sampling_freq 
"""
def phase_wavelet_cycles_min_ch(segment, frequencies, num_cycles, sampling_freq):
    
    T = segment.size
    phase_data = np.zeros((T, len(frequencies)))
    conv_data = np.zeros((T, len(frequencies)))
    wave_list, t_wave = [], []
    
    for f_idx, frequency in enumerate(frequencies):
        
        # Compute the wavelet 
        sigma = num_cycles / (2 * np.pi * frequency)   # s = n/(2*pi*f) -- standard deviation of the Gaussian
        t = np.arange(-3.5*sigma, 3.5*sigma, 1/sampling_freq)    # length of the wavelet 
        wavelet = np.exp(2j * np.pi * frequency * t) * np.exp(-t**2 / (2 * sigma**2))   # Morlet Wavelet 
        
        # store wavelets for each frequency
        wave_list.append(wavelet) 
        t_wave.append(t)
        
        # Convolve signal with wavelet
        convolution = convolve(segment, wavelet, 'same')
        # Extract instantaneous phase
        phase_data[:, f_idx] = np.angle(convolution)  # instantaneous phase, for each frequency, for each channel 
        conv_data[:, f_idx] = convolution.real
            
    return phase_data, wave_list, conv_data, t_wave


"""
Compute the instantaneous phase for the LFP chunks of data which are low speed.
- mask array contains 1 for low speed and 0 for high speed data. The LFP is masked with this
array for each channel and minute of the data, to get LFP low speed only.

LFP Low speed segments are identified and the instantaneous phase of these segments is computed 
separately (to avoid contributions in the convolution from high speed LFP, due to kernel outside
the low-speed LFP boundaries). These resulting segments of instantaneous phase are then stitched 
together by adding NaN in place of high speed segments.

The result is an multi-dim array (min, T, ch) with phase value for low speed LFP and NaN for high speed LFP


Inputs: 
lfp: lfp data dim: (min, T, ch) for the HPC pyramidal cells only
mask: low/high speed mask dim: (min, T, ch)
frequencies: range of frequencies for phase extraction (1D array)
n_cycle: n parameter in the Gaussian kernel for the number of cycles
fs = sampling frequency
"""

def phase_extraction_lfp_low_speed(lfp, mask, freq, n_cycles, fs):
    
    # Shape of the LFP signal
    min_dim, T_dim, ch_dim = lfp.shape
    
    # Placeholder for the final phase data, dim: (min, T, ch, freq)
    phase_data_full = np.full(lfp.shape + (len(freq),), np.nan) # phase 
    conv_data_full = np.full(lfp.shape + (len(freq),), np.nan)  # result of convolution (envelope at a given frequency)

    # Process each minute and channel separately
    for min_idx in range(min_dim):
        for ch_idx in range(ch_dim):
            
            # Extract the LFP signal for the current minute and channel
            current_lfp = lfp[min_idx, :, ch_idx]
            # Extract the corresponding mask with NaNs
            current_mask = mask[min_idx, :, ch_idx]

            # Find indices where the valid segments start and end, start will be 1, end will be -1
            change_points = np.diff(current_mask, prepend=0, append=0)
            start_indices = np.where(change_points == 1)[0]
            end_indices = np.where(change_points == -1)[0]

            # Process each valid segment
            for start, end in zip(start_indices, end_indices):
                
                # Extract the segment
                segment = current_lfp[start:end]
                
                # Compute the phase for the segment
                phase_segment, wave_list_segment, conv_segment, t_wave_segment = phase_wavelet_cycles_min_ch(segment, freq, n_cycles, fs)

                # Place the computed phase data in the corresponding location in the full phase data array
                phase_data_full[min_idx, start:end, ch_idx, :] = phase_segment    # (min, time, ch, freq)
                conv_data_full[min_idx, start:end, ch_idx, :] = conv_segment

    return phase_data_full, conv_data_full


"""
Compute the 2D map of frequency vs phase for the LFP-spike coupling 
Input: 
- phase of a single epoch (20 min long)
- mask for the same epoch for low speed LFP (20 min long)
- spike activity for the same epoch (20 min long)
Output:
- 2D map (freq, phase), with phase in [0, 360] degree
"""
def freq_phase_map(phase, mask_R, spk_epoch, ch, epoch_name):
    
    not_nan = np.where(mask_R[:,ch])[0] # Indexes for LFP low speed values
    n_freq = phase.shape[2] # number of frequencies 

    hist_list = []
    for f in range(n_freq):

        # phase values for not NaN (i.e. masked lfp) -- i.e. phase for low speed trials 
        phase_deg = np.degrees(phase[not_nan,ch,f])
        # Convert radiants to degrees in the range of 0 to 360
        phase_deg_adj = np.where(phase_deg < 0, phase_deg + 360, phase_deg)

        hist, bin_edges = np.histogram(phase_deg_adj, bins=36, range=[0, 360], weights=spk_epoch[not_nan])
        hist_list.append(hist)
        hist_arr = np.array(hist_list)
    print('spike count, low speed only {} {}'.format(epoch_name,np.sum(spk_epoch[not_nan])))
    return hist_arr, bin_edges

"""
Plot the frequency-phase map for one epoch only, i.e. baseline, low dose, etc..

Input: histogram_array matrix: each row contains the histogram of firing for 
each frequency
bin_edges: 1D array of bin values for the phase variable

Other parameters for plotting:
step: distance between phase values on x axis
step_f: distance between frequency values on y axis

"""

def plot_freq_phase_map_one_epoch(hist_arr, bin_edges, freq, title_name, step = 5, step_f = 4):


    row_totals = np.sum(hist_arr, axis=1) # normalization factor
    
    pdf = hist_arr / row_totals[:, np.newaxis]
    pdf_sm = gaussian_filter(np.nan_to_num(pdf), sigma=0.9) # smoothing 
    
    # Plotting the heatmap
    plt.figure(figsize=(3,5))
    plt.imshow(pdf_sm, cmap='RdYlBu_r', interpolation='nearest', origin='lower',aspect='auto')
    
    # xticks 
    bin_values = bin_edges[:-1]
    bins = bin_values[::step]
    xticks = np.arange(0,bin_edges[:-1].size,step)
    
    plt.xticks(ticks=xticks, labels=['{:d}'.format(int(b)) for b in bins ],rotation=45, fontsize=14)
    plt.yticks(ticks=np.arange(0, len(freq),step_f), labels=['{:d}'.format(int(f)) for f in freq[::step_f]],fontsize=14)
    cbar = plt.colorbar(label='Spike count')
    cbar.set_label('Probability', fontsize=12)  # Set a smaller font size for the label
    
    cbar.ax.tick_params(labelsize=12)
    
    # Adding axis labels with a custom font
    font = {'family': 'serif',
            'color':  'darkred',
            'weight': 'normal',
            'size': 12,
           }
    plt.xlabel('Phase (degrees)', fontdict=font)
    plt.ylabel('Frequency', fontdict=font)
    plt.title(title_name,fontsize=15)
    
    plt.show()



""" -------------------------------------------------------
Calculate normalized spiking count for each epoch
and Gaussian smoothed spiking count in the phase-frequency space
"""
# Define a function to calculate PDF and smoothed PDF
def calculate_pdf_sm(hist_arr, norm, sigma):
    pdf = hist_arr / norm # normalize histogram count with max count across epochs 
    pdf_sm = gaussian_filter(np.nan_to_num(pdf), sigma=sigma)
    return pdf, pdf_sm


""" -------------------------------------------------------
Plot subplots for the phase-frequency LFP-spike coupling for all the 4 epochs together, one subplot for each layer
"""
def setup_subplot(fig, ax, pdf, pdf_sm, freq, smooth, bin_edges, step_x, step_f, title, show_y_axis, bar_name, vmin, vmax):
    
    # Define a common font for titles
    title_font = {'family': 'sans-serif','fontsize': 14}
    label_font = {'family': 'sans-serif', 'color':  'black', 'weight': 'ultralight', 'size': 14}
    
    if smooth: distribution = pdf_sm 
    else: distribution = pdf
    
    # vmin = np.nanmin(distribution)
    # vmax = np.nanmax(distribution)
    
    # print('vmin local', vmin)
    # print('vmax local',vmax)
    
    cax = ax.imshow(1.*distribution, cmap='RdYlBu_r', interpolation='nearest', origin='lower', aspect='auto', vmin=vmin, vmax=vmax)
    bin_values = bin_edges[:-1]
    bins = bin_values[::step_x]
    
    xticks = np.arange(0, bin_edges[:-1].size, step_x)
    ax.set_xticks(xticks)
    ax.set_xticklabels(['' if i % 2!=0 else '{:d}'.format(int(b)) for i, b in enumerate(bins)], rotation=45, fontsize=14)
    ax.set_xlabel('Phase (deg)', fontdict=label_font)

    
    ax.set_yticks(np.arange(0, len(freq), step_f))
    if show_y_axis:
        ax.set_yticklabels(['{:d}'.format(int(f)) for f in freq[::step_f]], fontsize=14)
        ax.set_ylabel('Frequency (Hz)', fontdict=label_font)
    else:
        ax.set_yticklabels([])
        
    ax.spines['left'].set_linewidth(0.5)  # Adjust the left spine thickness
    ax.spines['bottom'].set_linewidth(0.5)  # Adjust the bottom spine thickness

    ax.tick_params(axis='both', which='major', width=1, length=7)  # Change '2' to your desired thickness
    
    

    # cbar = fig.colorbar(plt.cm.ScalarMappable(cmap='RdYlBu_r', norm=plt.Normalize(vmin=vmin, vmax=vmax)), ax=axes.ravel().tolist(), pad=0.01, aspect=20)
    # cbar = fig.colorbar(plt.cm.ScalarMappable(cmap='RdYlBu_r', norm=plt.Normalize(vmin=vmin, vmax=vmax)), pad=0.01, aspect=20)
    
    if bar_name: # print colorbar 
        cbar = fig.colorbar(cax, ax=ax, cmap='RdYlBu_r', format='%.3f')
        cbar.set_label('Probability', fontsize=12)
        cbar.ax.tick_params(labelsize=9,  width=1, length=7) 
        cbar.outline.set_linewidth(0.5)
        
        
""" -------------------------------------------------------
Plot frequency-phase plots for all the epochs together, one subplot for each layer

"""

def plot_freq_phase_map_all_epochs(hist_dict, bin_dict, norm_dict, freq, rec, sess, cell, idx, layer, layer_dir, layer_acr, lfp_ch, smooth = True, sigma = 1, save_flag = False, step_x = 5, step_f = 8):

    
    # norm_tot = norm_B + norm_L + norm_M + norm_H
    # Calculate PDFs
    pdf_B, pdf_sm_B = calculate_pdf_sm(hist_dict['B'], norm_dict['B'], sigma)
    pdf_L, pdf_sm_L = calculate_pdf_sm(hist_dict['L'], norm_dict['L'], sigma)
    pdf_M, pdf_sm_M = calculate_pdf_sm(hist_dict['M'], norm_dict['M'], sigma)
    pdf_H, pdf_sm_H = calculate_pdf_sm(hist_dict['H'], norm_dict['H'], sigma)
    
    # # Determine the global min and max values for the color scale
    
    vmin = min(np.nanmin(pdf_sm_B), np.nanmin(pdf_sm_L), np.nanmin(pdf_sm_M), np.nanmin(pdf_sm_H))
    vmax = max(np.nanmax(pdf_sm_B), np.nanmax(pdf_sm_L), np.nanmax(pdf_sm_M), np.nanmax(pdf_sm_H))
    print('vmin ', vmin, 'vmax', vmax)
    
    # Create 4 subplots
    fig, axes = plt.subplots(1, 4, figsize=(10, 5), constrained_layout=True) # Adjust the figsize as needed
    
    # Define a common font for titles
    title_font = {'family': 'sans-serif','fontsize': 14}
    label_font = {'family': 'sans-serif', 'color':  'black', 'weight': 'ultralight', 'size': 14}
    
    # Plot each histogram
    setup_subplot(fig, axes[0], pdf_B, pdf_sm_B, freq, smooth, bin_dict['B'], step_x, step_f, 'Baseline', True, False, vmin, vmax)
    setup_subplot(fig, axes[1], pdf_L, pdf_sm_L, freq, smooth, bin_dict['L'], step_x, step_f, 'Low dose', False, False, vmin, vmax)
    setup_subplot(fig, axes[2], pdf_M, pdf_sm_M, freq, smooth, bin_dict['M'], step_x, step_f, 'Mid dose', False, False, vmin, vmax)
    setup_subplot(fig, axes[3], pdf_H, pdf_sm_H, freq, smooth, bin_dict['H'], step_x, step_f, 'High dose',False, True, vmin, vmax)
    
    # Add a common color bar for all subplots
    # cbar = fig.colorbar(plt.cm.ScalarMappable(cmap='RdYlBu_r', norm=plt.Normalize(vmin=vmin, vmax=vmax)), ax=axes.ravel().tolist(), pad=0.01, aspect=20)
    # cbar.set_label('Probability', fontsize=12)
    # cbar.ax.tick_params(labelsize=10,  width=1, length=7) 
    # cbar.outline.set_linewidth(0.5)
    
    fig.suptitle(f'{layer}, RS Ket, cell = {cell}, lfp ch = {lfp_ch} ', fontsize=16, fontweight='regular')
    plt.show()
    
    if save_flag:
        save_figures_freq_phase(fig, rec, sess, 'HPC', cell, idx, layer_dir, layer_acr, lfp_ch)
        
    
""" -------------------------------------------------------
Save phase-frequency plot for the 4 epochs together 
"""
def save_figures_freq_phase(fig, rec, sess, brain_reg, cell, idx, layer_dir, layer_acr, lfp_ch):

    main_dir = r'C:\Users\fentonlab\Desktop\Gino\LFPs\HPC'
    path = rec[sess][brain_reg]

    dir_sess = path.split('\\')[-3]     # path for session directory
    full_dir_path = os.path.join(main_dir, dir_sess, f'Figures\\freq_phase\\all_channels\\cell_{idx}_ID_{cell}\\{layer_dir}')
    
    if not os.path.exists(full_dir_path):
        os.makedirs(full_dir_path)
    
    file_name = os.path.join(full_dir_path,'sess_{}_{}_cell_{}_lfp_ch_{}_freq_phase_300_Hz_CSD.pdf'.format(sess, layer_acr, cell, lfp_ch))
    fig.savefig(file_name, dpi=300)
    file_name = os.path.join(full_dir_path,'sess_{}_{}_cell_{}_lfp_ch_{}_freq_phase_300_Hz_CSD.png'.format(sess, layer_acr, cell, lfp_ch))
    fig.savefig(file_name, dpi=300)
    print("\nsaving file:\n",file_name)
    
    
""" -------------------------------------------------------
Save Spike-Count histogram
"""
def save_spike_count_hist(fig, rec, sess, brain_reg, cell, layer_dir, layer_acr):

    main_dir = r'C:\Users\fentonlab\Desktop\Gino\LFPs\HPC'
    path = rec[sess][brain_reg]

    dir_sess = path.split('\\')[-3]     # path for session directory
    full_dir_path = os.path.join(main_dir, dir_sess, 'Figures\\freq_phase\\all_channels', layer_dir)
    
    if not os.path.exists(full_dir_path):
        os.makedirs(full_dir_path)
    
    file_name = os.path.join(full_dir_path,'sess_{}_{}_cell_{}_lfp_ch_{}_spike_count_hist_CSD.pdf'.format(sess, layer_acr, cell, lfp_ch))
    print(file_name)
    fig.savefig(file_name, dpi=300)
    file_name = os.path.join(full_dir_path,'sess_{}_{}_cell_{}_lfp_ch_{}_spike_count_hist_CSD.png'.format(sess, layer_acr, cell, lfp_ch))
    fig.savefig(file_name, dpi=300)
    
    
        
""" -------------------------------------------------------
Plot subplots for the phase-frequency LFP-spike coupling for all the 4 epochs together
"""

def setup_subplot_4_by_4(fig, ax, pdf, pdf_sm, freq, smooth, bin_edges, step_x, step_f, title, i, title_print, stratus_name, show_y_axis, bar_name, vmin, vmax):
    
    # Define a common font for titles
    title_font = {'family': 'Arial','fontsize': 14}
    label_font = {'family': 'Arial', 'color':  'black', 'weight': 'ultralight', 'size': 14}
    
    if smooth: distribution = pdf_sm 
    else: distribution = pdf
    
    # vmin = np.nanmin(distribution)
    # vmax = np.nanmax(distribution)
    
    # print('vmin local', vmin)
    # print('vmax local',vmax)
    
    cax = ax.imshow(1.*distribution, cmap='RdYlBu_r', interpolation='nearest', origin='lower', aspect='auto', vmin=vmin, vmax=vmax)
    bin_values = bin_edges[:-1]
    bins = bin_values[::step_x]
    
    xticks = np.arange(0, bin_edges[:-1].size, step_x)
    ax.set_xticks(xticks)
    if i == 3: # if you are plotting the last layer at the bottom, add x labels 
        ax.set_xticklabels(['' if i % 2!=0 else '{:d}'.format(int(b)) for i, b in enumerate(bins)], rotation=45, fontsize=14)
        ax.set_xlabel('Phase (deg)', fontdict=label_font)
    else:
        ax.set_xticklabels([])
    
    ax.set_yticks(np.arange(0, len(freq), step_f))
    if show_y_axis:
        ax.set_yticklabels(['{:d}'.format(int(f)) for f in freq[::step_f]], fontsize=14)
        ax.set_ylabel(f'{stratus_name.upper()} \n Frequency (Hz)', fontdict=label_font)
    else:
        ax.set_yticklabels([])
        
    ax.spines['left'].set_linewidth(0.5)  # Adjust the left spine thickness
    ax.spines['bottom'].set_linewidth(0.5)  # Adjust the bottom spine thickness

    ax.tick_params(axis='both', which='major', width=1, length=7)  # Change '2' to your desired thickness
    
    if title_print == 0: # set suplots titles, i.e. baseline, low, mid, high dose
        ax.set_title(title, fontdict=label_font)
    
    # cbar = fig.colorbar(plt.cm.ScalarMappable(cmap='RdYlBu_r', norm=plt.Normalize(vmin=vmin, vmax=vmax)), ax=axes.ravel().tolist(), pad=0.01, aspect=20)
    # cbar = fig.colorbar(plt.cm.ScalarMappable(cmap='RdYlBu_r', norm=plt.Normalize(vmin=vmin, vmax=vmax)), pad=0.01, aspect=20)
    
    if bar_name: # print colorbar 
        cbar = fig.colorbar(cax, ax=ax, cmap='RdYlBu_r', format='%.3f')
        cbar.set_label('Probability', fontsize=14, fontdict=label_font)
        cbar.ax.tick_params(labelsize=12,  width=1, length=7) 
        cbar.outline.set_linewidth(0.5)
    
""" -------------------------------------------------------
Plot spiking histogram (counts) across epochs -- Not in use for now, 1 plot for each epoch.
"""
def plot_save_spiking_count_across_epochs(norm_dict, idx, cell, stratum_name):
        
        
    spike_count = [norm_dict['B'], norm_dict['L'], norm_dict['M'], norm_dict['H']] # spiking count for each epoch
    bins = range(len(spike_count))
    
    # Plot bars
    fig = plt.figure(figsize=(4,4))
    plt.bar(bins, spike_count, edgecolor='white', color='lightgrey', align='center')
    
    plt.title('Spiking count across epochs, {}, cell: {} ID: {}'.format(stratum_name.upper(),idx, cell), fontsize=10)
    plt.xlabel('Epoch',fontsize=14)
    plt.ylabel('Count',fontsize=14)
    plt.xticks(ticks=bins, labels=["base","low","mid","high"],fontsize=12)
    plt.yticks(fontsize=12)
    
    for axis in ['top','bottom','left','right']:
        plt.gca().spines[axis].set_linewidth(0.5)
    plt.gca().tick_params(axis='x', length=4, width=1)
    plt.gca().tick_params(axis='y', length=4, width=1)
    
    plt.tight_layout()
    plt.show()
    
    # save_spike_count_hist(fig, rec, sess, 'HPC', cell, layer_dir, layer_acr)
    
        

""" -------------------------------------------------------
Plot frequency-phase plots for all the strata together, in a 4x4 grid plot, where on the x we have: baseline, low, mid, high dose
on the y we have stratum oriens, pyramidale, radiatum, locmol 
"""


def plot_save_freq_phase_map_all_epochs_4_by_4(hist_dict_array, bin_dict_array, norm_dict_array, freq, rec, sess, cell, idx, stratum_name, 
                                          smooth=True, sigma=1, save_flag=False, step_x=5, step_f=8):
    title_font = {'family': 'Arial','size': 14}
    label_font = {'family': 'Arial', 'color':  'black', 'weight': 'ultralight', 'size': 14}
    pdf_sm_dict = [] # list of dictionaries to store the phase-freq map for each epoch, and each stratum 
    
    
    fig, axes = plt.subplots(4, 4, figsize=(10, 16), constrained_layout=True)  # Create a 4x4 grid of subplots
    
    # Plot each layer on a single row, made of 4 columns 
    for ch, (hist_dict, bin_dict, norm_dict) in enumerate(zip(hist_dict_array, bin_dict_array, norm_dict_array)):
        # Calculate row-specific vmin and vmax for color scaling
        vmin, vmax = calculate_row_vmin_vmax(hist_dict, norm_dict, sigma)
                   
        
        # Plot each histogram
        pdf_B, pdf_sm_B = calculate_pdf_sm(hist_dict['B'], norm_dict['B'], sigma)
        setup_subplot_4_by_4(fig, axes[ch,0], pdf_B, pdf_sm_B, freq, smooth, bin_dict['B'], step_x, step_f, 'Baseline', ch, ch, stratum_name[ch], True, False, vmin, vmax)
        
        pdf_L, pdf_sm_L = calculate_pdf_sm(hist_dict['L'], norm_dict['L'], sigma)
        setup_subplot_4_by_4(fig, axes[ch,1], pdf_L, pdf_sm_L, freq, smooth, bin_dict['L'], step_x, step_f, 'Low dose', ch, ch, stratum_name[ch], False, False,vmin, vmax)
        
        pdf_M, pdf_sm_M = calculate_pdf_sm(hist_dict['M'], norm_dict['M'], sigma)
        setup_subplot_4_by_4(fig, axes[ch,2], pdf_M, pdf_sm_M, freq, smooth, bin_dict['M'], step_x, step_f, 'Mid dose', ch, ch, stratum_name[ch], False, False,vmin, vmax)
        
        pdf_H, pdf_sm_H = calculate_pdf_sm(hist_dict['H'], norm_dict['H'], sigma)
        setup_subplot_4_by_4(fig, axes[ch,3], pdf_H, pdf_sm_H, freq, smooth, bin_dict['H'], step_x, step_f, 'High dose', ch, ch, stratum_name[ch], False, True, vmin, vmax)
    
        pdf_sm_dict.append({'B':pdf_sm_B, 'L':pdf_sm_L, 'M':pdf_sm_M, 'H':pdf_sm_H}) # append freq-maps for eahc stratum (all epochs)
        
        
        
    fig.suptitle(f'RS Ket, cell: {idx}, id: {cell}, CA1 layers ', fontweight='regular',fontdict=title_font)
    plt.show()
    
    if save_flag:
        save_figures_freq_phase_4_by_4(fig, rec, sess, 'HPC', cell, idx) # save figures for each cell
        save_matrices_freq_phase_4_by_4(pdf_sm_dict, rec, sess, 'HPC', cell, idx) # save matrices freq-phase for each cell 
        pass



""" -------------------------------------------------------
Save phase-frequency plot for the 4 epochs together 
"""
def save_figures_freq_phase_4_by_4(fig, rec, sess, brain_reg, cell, idx):

    main_dir = r'C:\Users\fentonlab\Desktop\Gino\LFPs\HPC'
    path = rec[sess][brain_reg]

    dir_sess = path.split('\\')[-3]     # path for session directory
    full_dir_path = os.path.join(main_dir, dir_sess, 'Figures\\freq_phase\\all_strata')
    
    if not os.path.exists(full_dir_path):
        os.makedirs(full_dir_path)
    
    file_name = os.path.join(full_dir_path,'sess_{}_cell_{}_id_{}_freq_phase_CSD.pdf'.format(sess,idx, cell))
    fig.savefig(file_name, dpi=300)
    file_name = os.path.join(full_dir_path,'sess_{}_cell_{}_id_{}_freq_phase_CSD.png'.format(sess,idx, cell))
    fig.savefig(file_name, dpi=300)
    print("\nsaving Figures on file:\n",file_name)
    
    
    
    
""" -------------------------------------------------------
Save phase-frequency matrices for the 4 epochs together 
in a dictionary
output: dictionary saved on file with the pdf_smoothed for a given cell-LFP
"""
def save_matrices_freq_phase_4_by_4(pdf_sm_dict, rec, sess, brain_reg, cell, idx):
    
    import pickle 
    
    main_dir = r'C:\Users\fentonlab\Desktop\Gino\LFPs\HPC'
    path = rec[sess][brain_reg]

    dir_sess = path.split('\\')[-3]     # path for session directory
    full_dir_path = os.path.join(main_dir, dir_sess, 'freq_phase_matrices')
    
    if not os.path.exists(full_dir_path):
        os.makedirs(full_dir_path)
    
    file_name = os.path.join(full_dir_path,'sess_{}_cell_{}_id_{}_freq_phase_CSD.pkl'.format(sess,idx, cell))
    with open(file_name,'wb')as file:
        pickle.dump(pdf_sm_dict, file)

    print("\nsaving matrices on file:\n",file_name)


""" -------------------------------------------------------
Load phase-frequency matrices for the 4 epochs together 
in a dictionary
"""
    
def load_matrices_freq_phase_4_by_4(rec, sess, brain_reg, cell, idx):
    
    import pickle 
    
    main_dir = r'C:\Users\fentonlab\Desktop\Gino\LFPs\HPC'
    path = rec[sess][brain_reg]

    dir_sess = path.split('\\')[-3]     # path for session directory
    full_dir_path = os.path.join(main_dir, dir_sess, 'freq_phase_matrices')
    
    if not os.path.exists(full_dir_path):
        os.makedirs(full_dir_path)
    
    file_name = os.path.join(full_dir_path,'sess_{}_cell_{}_id_{}_freq_phase_CSD.pkl'.format(sess,idx, cell))
    with open(file_name,'rb')as file:
        pdf_sm_dict = pickle.load(file)

    return pdf_sm_dict
    
""" -------------------------------------------------------
Plot spiking histogram (counts) across epochs and strata - One figure, 4 plots 
"""

def plot_save_spiking_count_across_epochs_and_strata(norm_dict_arr, rec, sess, idx, cell, stratum_name):
    
    fig, axs = plt.subplots(4, 1, figsize=(4,16))  # 4 rows, 1 column
    
    for i in range(4): # for each stratum
    
        spike_count = [norm_dict_arr[i]['B'], norm_dict_arr[i]['L'], norm_dict_arr[i]['M'], norm_dict_arr[i]['H']]
        bins = range(len(spike_count))
        
        axs[i].bar(bins, spike_count, edgecolor='white', color='lightgrey', align='center')
        
        axs[i].set_title('Spiking count, {}, cell: {} ID: {}'.format(stratum_name[i].upper(),idx, cell), fontsize=10)
        axs[i].set_xlabel('Epoch', fontsize=14)
        axs[i].set_ylabel('Count', fontsize=14)
        axs[i].set_xticks(ticks=bins)
        axs[i].set_xticklabels(labels=["base","low","mid","high"], fontsize=12)
        axs[i].tick_params(axis='x', length=4, width=1, labelsize=12)
        axs[i].tick_params(axis='y', length=4, width=1, labelsize=12)
        
        for axis in ['top','bottom','left','right']:
            axs[i].spines[axis].set_linewidth(0.5)
    
    plt.tight_layout()
    plt.show()
    
    save_figures_spiking_histo(fig, rec, sess, "HPC", cell, idx)
    

""" -------------------------------------------------------
Save spiking count plot for all epochs, all strata 
"""
def save_figures_spiking_histo(fig, rec, sess, brain_reg, cell, idx):

    main_dir = r'C:\Users\fentonlab\Desktop\Gino\LFPs\HPC'
    path = rec[sess][brain_reg]

    dir_sess = path.split('\\')[-3]     # path for session directory
    full_dir_path = os.path.join(main_dir, dir_sess, 'Figures\\freq_phase\\all_strata\\spiking_histo')
    
    if not os.path.exists(full_dir_path):
        os.makedirs(full_dir_path)
    
    file_name = os.path.join(full_dir_path,'sess_{}_cell_{}_id_{}_spiking_histo.pdf'.format(sess,idx, cell))
    fig.savefig(file_name, dpi=300)
    file_name = os.path.join(full_dir_path,'sess_{}_cell_{}_id_{}_spiking_histo.png'.format(sess,idx, cell))
    fig.savefig(file_name, dpi=300)
    print("\nsaving Figures on file:\n",file_name)
    

""" -------------------------------------------------------
Calculate max and min in frequency-phase histogram 
"""

def calculate_row_vmin_vmax(hist_dict, norm_dict, sigma):
    vmin = vmax = None
    for key in ['B', 'L', 'M', 'H']:
        pdf, pdf_sm = calculate_pdf_sm(hist_dict[key], norm_dict[key], sigma)
        if vmin is None or np.nanmin(pdf_sm) < vmin:
            vmin = np.nanmin(pdf_sm)
        if vmax is None or np.nanmax(pdf_sm) > vmax:
            vmax = np.nanmax(pdf_sm)
    return vmin, vmax


""" -------------------------------------------------------
Average LFP across stratum 
"""

def get_lfp_stratum(lfp_B, lfp_L, lfp_M, lfp_H, stratum_chs):
    
    # get the average LFP for each stratum 
    lfp_B_str = np.mean(lfp_B[:,:,stratum_chs], axis=2) # baseline
    lfp_L_str = np.mean(lfp_L[:,:,stratum_chs], axis=2)  # low dosage
    lfp_M_str = np.mean(lfp_M[:,:,stratum_chs], axis=2) # mid dosage
    lfp_H_str = np.mean(lfp_H[:,:,stratum_chs], axis=2) # high dosage    
    
    return lfp_B_str, lfp_L_str, lfp_M_str, lfp_H_str



""" -------------------------------------------------------
Combine masks across channels for the same CA1 stratum 
"""

def get_mask_stratum(mask_B, mask_L, mask_M, mask_H, stratum_chs):

    min_dim = mask_B.shape[0]
    mask_B_list, mask_L_list, mask_M_list, mask_H_list = [], [], [], []
    
    for minute in range(min_dim):
        
        mask_minute = np.all(mask_B[minute,:,stratum_chs],axis=0).astype(int)
        mask_B_list.append(mask_minute)
        mask_minute = np.all(mask_L[minute,:,stratum_chs],axis=0).astype(int)
        mask_L_list.append(mask_minute)
        mask_minute = np.all(mask_M[minute,:,stratum_chs],axis=0).astype(int)
        mask_M_list.append(mask_minute)
        mask_minute = np.all(mask_H[minute,:,stratum_chs],axis=0).astype(int)
        mask_H_list.append(mask_minute)
        
    return np.array(mask_B_list), np.array(mask_L_list), np.array(mask_M_list), np.array(mask_H_list)



def mean_abs_error_freq_phase_maps(pdf_sm,theta,low_gamma,high_gamma):
    
    MAE = []
    for ch in range(4):
        MAE.append({
            'theta': {
                'BL': np.mean(np.abs(pdf_sm[ch]["B"][theta] - pdf_sm[ch]["L"][theta])), # baseline-low dose
                'BM': np.mean(np.abs(pdf_sm[ch]["B"][theta] - pdf_sm[ch]["M"][theta])), # baseline-mid dose
                'BH': np.mean(np.abs(pdf_sm[ch]["B"][theta] - pdf_sm[ch]["H"][theta])), # baseline-high dose
                }, 
            'low_gamma': {
                'BL': np.mean(np.abs(pdf_sm[ch]["B"][low_gamma] - pdf_sm[ch]["L"][low_gamma])),
                'BM': np.mean(np.abs(pdf_sm[ch]["B"][low_gamma] - pdf_sm[ch]["M"][low_gamma])),
                'BH': np.mean(np.abs(pdf_sm[ch]["B"][low_gamma] - pdf_sm[ch]["H"][low_gamma])),
                },
            'high_gamma': {
                'BL': np.mean(np.abs(pdf_sm[ch]["B"][high_gamma] - pdf_sm[ch]["L"][high_gamma])),
                'BM': np.mean(np.abs(pdf_sm[ch]["B"][high_gamma] - pdf_sm[ch]["M"][high_gamma])),
                'BH': np.mean(np.abs(pdf_sm[ch]["B"][high_gamma] - pdf_sm[ch]["H"][high_gamma])),
                }
            })
        
        
    return  MAE
    
""""
Prepare data to 3D plot, each axis is baseline-low, baseline-mid, baseline-high
ch: channel of the stratum, i.e. 4 channels in total
freq_band: either theta, low_gamma, high_gamma
"""
def create_array_for_3D_plots(mae_sess,idx_cell_HPC,ch,freq_band):
    
    vector = []
    for idx in range(len(idx_cell_HPC)):
        
        x = mae_sess[idx][ch][freq_band]['BL']
        y = mae_sess[idx][ch][freq_band]['BM']
        z = mae_sess[idx][ch][freq_band]['BH']
        
        vector.append([idx,x,y,z])
    
    
    return np.array(vector)


"""
create array for 3D plot, with axis theta, low_gamma, high)gamma.
Output: return a 3D vector with the mean absolute value of the epoch_diff
theta, low gamma, high gamma
"""

def create_array_for_3D_plots_axis_freq_band(mae_sess,idx_cell_HPC,ch,epoch_diff):
    
    vector = []
    for idx in range(len(idx_cell_HPC)):
        
        theta = mae_sess[idx][ch]['theta'][epoch_diff]
        low_gamma = mae_sess[idx][ch]["low_gamma"][epoch_diff]
        high_gamma = mae_sess[idx][ch]["high_gamma"][epoch_diff]
        
        vector.append([idx,theta,low_gamma,high_gamma])
    
    
    return np.array(vector)


"""
3D plot, with axis theta, low_gamma, high)gamma
    for each cell in the session. Plot the mean absolute error 
    between baseline-ketamine dose
"""
    
def plot_3D_mean_abs_error_all_cells(vector, title_name):

    X = vector[:,1]
    Y = vector[:,2]
    Z = vector[:,3]
    
    d = 0.005
    
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')
    
    # indices = np.arange(vector.shape[0])
    # norm = (indices - indices.min()) / (indices.max() - indices.min())
    # # Choose a colormap
    # cmap = plt.get_cmap('viridis')
    # colors = cmap(norm)

    # Scatter plot with improved appearance
    # Color points with a gradient based on the Z value, change size with s parameter
    scatter = ax.scatter(X, Y, Z, c='k', s=50, edgecolor='k', alpha=0.7)
    
    for i, (x, y, z) in enumerate(zip(X, Y, Z)):
        ax.text(x, y, z+d, f'{i}', color='black',fontsize=12)
    
    # Enhancing labels with font size adjustments
    ax.set_xlabel('theta', fontsize=15, fontweight='light')
    ax.set_ylabel('low-gamma', fontsize=15, fontweight='light')
    ax.set_zlabel('high-gamma', fontsize=15, fontweight='light')
    
    # Setting the title with additional styling
    ax.set_title(title_name, fontsize=18, fontweight='light', color='k')
    
    # Setting the background color of the plot
    ax.set_facecolor('white')  # Change background color
    
    # Configuring the grid
    ax.grid(True, linestyle='--', color='gray', linewidth=0.3, alpha=0.8)
    # Make the numbers on the axes smaller
    ax.tick_params(axis='both', which='major', labelsize=10)
    ax.tick_params(axis='both', which='minor', labelsize=6)
        
    # Show plot with improved appearance
    plt.show()
    

        
def plot_3D_mean_abs_error_all_cells_animated(vector, title_name):

    from matplotlib.animation import FuncAnimation
    from mpl_toolkits.mplot3d import Axes3D
    
    X = vector[:,1]
    Y = vector[:,2]
    Z = vector[:,3]
    
    d = 0.003
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    scatter = ax.scatter(X, Y, Z, c='k', s=50, edgecolor='k', alpha=0.7)
    
    ax.set_xlabel('theta', fontsize=15, fontweight='light')
    ax.set_ylabel('low-gamma', fontsize=15, fontweight='light')
    ax.set_zlabel('high-gamma', fontsize=15, fontweight='light')
    
    # Setting the title with additional styling
    ax.set_title(title_name, fontsize=18, fontweight='light', color='k')
    
    # Setting the background color of the plot
    ax.set_facecolor('white')  # Change background color
    
    # Configuring the grid
    ax.grid(True, linestyle='--', color='gray', linewidth=0.3, alpha=0.8)
    # Make the numbers on the axes smaller
    ax.tick_params(axis='both', which='major', labelsize=10)
    ax.tick_params(axis='both', which='minor', labelsize=6)
    
    # Setting the initial view
    ax.view_init(elev=20, azim=30)
    
    # Update function for the animation
    def update(frame):
        # Rotating by updating the azimuth angle (azim)
        ax.view_init(elev=20, azim=frame)
        return scatter,
    
    # Creating the animation
    ani = FuncAnimation(fig, update, frames=np.arange(0, 360, 2), blit=False)

    # To display the animation inline in a Jupyter notebook, use:
    from IPython.display import HTML
    HTML(ani.to_html5_video())
    # Save the animation
    # ani.save('3d_animation.mp4', writer='ffmpeg')
    
    # # Show the plot
    # plt.show()

    
    
    
    
    