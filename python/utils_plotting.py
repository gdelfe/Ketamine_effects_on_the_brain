# -*- coding: utf-8 -*-
"""
Spyder Editor

@ Gino Del Ferraro, Fenton lab, Oct 2023
"""

import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# PLOT SPEED FUNCTIONS 
# =============================================================================

def plot_speed(speed_up,sess,start,end,Xtick,fs_lfp=2500,fs_beh=100):

    speed_T_up = np.linspace(0,speed_up.size,speed_up.size) # t axis for speed
    
    X = Xtick   # interval for xticks
    s = fs_lfp*start # star
    e = fs_lfp*end # end
    
    print('Number of seconds: ',(e-s)/fs_lfp)
    
    plt.figure(figsize=(12,5))
    
    plt.plot(speed_T_up[s:e],speed_up[s:e],linewidth=0.5,color='orange') # upsampled signal 
    
    plt.title('Speed, given window',fontsize=12)
    plt.xlabel('time (sec)',fontsize=12)
    plt.ylabel('speed (mm/s)',fontsize=12)
    
    xticks = np.arange(0,e-s, step=X*fs_lfp)
    yticks = np.arange(0,400,50)
    plt.xticks(ticks=xticks, labels=[str(int(i/fs_lfp)) for i in xticks],fontsize=10)
    plt.yticks(ticks=yticks,fontsize=10)
    plt.grid(True,linewidth=0.2)
    plt.ylim([0,200])
             
    plt.tight_layout()
    plt.show()
    
    plt.figure(figsize=(12,5))
    plt.plot(speed_up[:],linewidth=0.5) # upsampled signal 
    plt.title('Speed whole recording',fontsize=12)
    plt.xlabel('time (min)',fontsize=12)
    plt.ylabel('speed (mm/s)',fontsize=12)
    xticks = np.arange(0,speed_up.size, step=10*60*2500) # every 10 min
    plt.xticks(ticks=xticks, labels=[str(int(i/60/2500)) for i in xticks],fontsize=10)
    plt.yticks(ticks=yticks,fontsize=10)
    plt.grid(True,linewidth=0.2)
    # plt.ylim([0,100])
    
    plt.tight_layout()
    plt.show()

# =============================================================================


def plot_speed_histo_logscale(speed,sess):
    
    plt.hist(speed[sess], bins=50, edgecolor = 'white', color= 'skyblue',density=True)
    plt.yscale('log')
    plt.xscale('linear')
    plt.title('speed histogram Log scale',fontsize=12)
    plt.xlabel('Value',fontsize=12)
    plt.ylabel('count',fontsize=12)
    plt.show()

# =============================================================================


def plot_speed_histo_regimes(speed_up, th_low = 50,th_mid = 100):
    

    low = np.sum(speed_up < th_low)
    mid = np.sum((speed_up >= th_low) & (speed_up <th_mid))
    high = np.sum(speed_up >=th_mid)
    speed_count = [low,mid,high]
    
    print('speed count ',speed_count)


    plt.figure(figsize=(4,4))
    bins = range(len(speed_count))
    
    plt.figure(figsize = (3,3))
    # Plot bars
    plt.bar(bins, speed_count, edgecolor='white', color='skyblue', align='center')
    
    plt.title('Speed histogram, 3 regimes',fontsize=12)
    plt.xlabel('speed regime',fontsize=12)
    plt.ylabel('Count',fontsize=12)
    plt.xticks(ticks=bins, labels=["low","mid","high"],fontsize=12)
    plt.yticks(fontsize=12)
    
    plt.show()



# =============================================================================
# PLOT LFP FUNCTIONS
# =============================================================================

def plot_lfp_two_channels(Lfp,ch1,ch2,start,end,Xtick,N=2500):

    L_start = start*N # start of time series in sec
    L_end = end*N # end of time series 
    T = Xtick*N # xtick period 
    
    
    plt.figure(figsize=(10,5))
    plt.plot(Lfp[L_start:L_end,ch1],linewidth=0.5)
    plt.title(f'LFP for channel {ch1}',fontsize=12)
    plt.xlabel('time (sec)',fontsize=12)
    plt.ylabel('Unknown UOM',fontsize=12)
    
    
    xticks = np.arange(0,L_end-L_start, step=T)
    plt.xticks(ticks=xticks, labels=[str(i/N) for i in xticks],fontsize=10)
    plt.yticks(fontsize=10)
    # print(xticks)
    plt.tight_layout()
    plt.show()
    
    
    plt.figure(figsize=(10,5))
    plt.plot(Lfp[L_start:L_end,ch2],linewidth=0.5)
    plt.title(f'LFP for channel {ch2}',fontsize=12)
    plt.xlabel('time (sec)',fontsize=12)
    plt.ylabel('Unknown UOM',fontsize=12)
    
    
    xticks = np.arange(0,L_end-L_start, step=T)
    plt.xticks(ticks=xticks, labels=[str(i/N) for i in xticks],fontsize=10)
    plt.yticks(fontsize=10)
    
    plt.tight_layout()
    plt.show()

# =============================================================================


def plot_lfp_two_channels_together(Lfp,ch1,ch2,start,end,Xtick,N=2500):

    L_start = start*N # start of time series in sec
    L_end = end*N # end of time series 
    T = Xtick*N # xtick period 
    
    
    plt.figure(figsize=(10,5))
    plt.plot(Lfp[L_start:L_end,ch1],linewidth=0.5)
    plt.plot(Lfp[L_start:L_end,ch2],linewidth=0.5)
    plt.title(f'LFP for channel {ch1} and bad channel {ch2}',fontsize=12)
    plt.xlabel('time (sec)',fontsize=12)
    plt.ylabel('Unknown UOM',fontsize=12)
    
    
    xticks = np.arange(0,L_end-L_start, step=T)
    plt.xticks(ticks=xticks, labels=[str(i/N) for i in xticks],fontsize=10)
    plt.yticks(fontsize=10)
    # print(xticks)
    plt.tight_layout()
    plt.show()
    


# =============================================================================


# P rows, Q columns 

def plot_lfp_various_channels(Lfp,ch1,ch2,start,end,P,Q,Xtick,N=2500):
    
    
    L_start = start*N # start of time series in sec
    L_end = end*N # end of time series 
    T = Xtick*N # xtick period 
    
    
    # uncomment for zscore spectrogram 
    # zlfp = zscore(Lfp[L_start:L_end,:],axis=0)
    # lfp_plot = zlfp
    
    lfp_plot = Lfp[L_start:L_end,:]
    
    plt.figure(figsize=(10,12))
    i = 1
    for ch in range(ch1,ch2):
        
        plt.subplot(P,Q,i)
        i += 1   
        
        plt.plot(lfp_plot[:,ch],linewidth=0.5)
        plt.title(f'LFP for channel {ch}',fontsize=12)
        plt.xlabel('time (sec)',fontsize=12)
        plt.ylabel('Unknown UOM',fontsize=12)
    
    
        xticks = np.arange(0,L_end-L_start, step=T)
        plt.xticks(ticks=xticks, labels=[str(i/N) for i in xticks],fontsize=10)
        plt.yticks(fontsize=10)
        plt.ylim([-100,100])
    
    plt.tight_layout()
    plt.show()

# =============================================================================


def plot_filtered_lfp(lfp_filt, start, end, Xtick, ch, N = 2500):
    
    
    L_start = start*N # start of time series in sec
    L_end = end*N # end of time series 
    T = Xtick*N # xtick period 
    
    
    plt.figure(figsize=(10,5))
    plt.plot(lfp_filt[L_start:L_end,ch],linewidth=0.5)
    # plt.plot(lfp_filt_B2[L_start:L_end,ch],linewidth=0.5)
    plt.title(f'LFP for channel {ch}',fontsize=12)
    plt.xlabel('time (sec)',fontsize=12)
    plt.ylabel('$\mu$V',fontsize=12)
    
    
    xticks = np.arange(0,L_end-L_start, step=T)
    plt.xticks(ticks=xticks, labels=[str(i/N) for i in xticks],fontsize=10)
    plt.yticks(fontsize=10)
    # print(xticks)
    plt.tight_layout()
    plt.show()
    
    