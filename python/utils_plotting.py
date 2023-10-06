# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 14:27:03 2023

@author: fentonlab
"""
import numpy as np
import matplotlib.pyplot as plt


def plot_speed(speed_up,sess,start,end,Xtick,fs_lfp=2500,fs_beh=100):

    speed_T_up = np.linspace(0,speed_up.size,speed_up.size) # t axis for speed
    
    X = fs_lfp*Xtick   # interval for xticks
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


def plot_speed_histo_logscale(speed,sess):
    
    plt.hist(speed[sess], bins=50, edgecolor = 'white', color= 'skyblue',density=True)
    plt.yscale('log')
    plt.xscale('linear')
    plt.title('speed histogram Log scale',fontsize=12)
    plt.xlabel('Value',fontsize=12)
    plt.ylabel('count',fontsize=12)
    plt.show()


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


