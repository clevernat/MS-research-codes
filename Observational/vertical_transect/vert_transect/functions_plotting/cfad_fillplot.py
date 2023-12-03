import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
from matplotlib.transforms import Bbox
from matplotlib.lines import Line2D

def cfad_fillplot(ref_s,ref_w,vel_s,vel_w,sw_s, sw_w,bin, height, h, id):

    ref_ave_w = np.zeros(h)
    ref_ave_s = np.zeros(h)
    vel_ave_w = np.zeros(h)
    vel_ave_s = np.zeros(h)
    sw_ave_w = np.zeros(h)
    sw_ave_s = np.zeros(h)
    ref_trans_w = np.transpose(ref_w)
    ref_trans_s = np.transpose(ref_s)
    vel_trans_w = np.transpose(vel_w)
    vel_trans_s = np.transpose(vel_s)
    sw_trans_w = np.transpose(sw_w)
    sw_trans_s = np.transpose(sw_s)

    for j in range(0, h):
        ref_ave_w[j] = np.nanmean(10**(ref_trans_w[j]/10.0))
        ref_ave_s[j] = np.nanmean(10**(ref_trans_s[j]/10.0))
        vel_ave_w[j] = np.nanmean(vel_trans_w[j])
        vel_ave_s[j] = np.nanmean(vel_trans_s[j])
        sw_ave_w[j] = np.nanmean(sw_trans_w[j])
        sw_ave_s[j] = np.nanmean(sw_trans_s[j])
        ref_ave_w[j] = 10.0*np.log10(ref_ave_w[j])
        ref_ave_s[j] = 10.0*np.log10(ref_ave_s[j])

    fig, ax = plt.subplots(1,3, figsize=(30,6))
    y = np.linspace(height[0],height[-1],len(height)+1)
    ax[0].plot(ref_ave_w, height, c = 'blue',linewidth=2.3)
    ax[0].plot(ref_ave_s, height, c = 'red',linewidth=2.3)
    ax[0].set_ylabel("MSL [km]",fontsize = 15)
    ax[0].set_xlabel("Reflectivity",fontsize = 15)
    ax[0].set_ylim(0,6000)
    ax[0].tick_params(axis='both', which='major', labelsize=15)
    ax[0].set_yticks([0,1000,2000,3000,4000,5000,6000])
    ax[0].set_yticklabels([0,1,2,3,4,5,6],fontsize = 15)
    ax[0].set_xlim(-30,30)
    ax[0].set_xticks(np.arange(-30,31,10))
    ax[0].set_xticklabels(np.arange(-30,31,10), fontsize =15)

    ax[1].plot(vel_ave_w, height, c = 'blue',linewidth=2.3)
    ax[1].plot(vel_ave_s, height, c = 'red',linewidth=2.3)
    ax[1].set_ylim(0,6000)
    ax[1].tick_params(axis='both', which='major', labelsize=15)
    ax[1].set_yticks([0,1000,2000,3000,4000,5000,6000])
    ax[1].set_xlabel("Doppler Velocity",fontsize = 15)
    ax[1].set_yticklabels([0,1,2,3,4,5,6],fontsize = 15)
    ax[1].set_xlim(-2.5,1.5)
    ax[1].set_xticks(np.round(np.arange(-2.5,1.6,0.5),1))
    ax[1].set_xticklabels(np.round(np.arange(-2.5,1.6,0.5),1), fontsize =15)

    ax[2].plot(sw_ave_w, height, c = 'blue',linewidth=2.3)
    ax[2].plot(sw_ave_s, height, c = 'red',linewidth=2.3)
    ax[2].set_ylim(0,6000)
    ax[2].tick_params(axis='both', which='major', labelsize=15)
    ax[2].set_yticks([0,1000,2000,3000,4000,5000,6000])
    ax[2].set_xlabel("Spectral Width",fontsize = 15)
    ax[2].set_yticklabels([0,1,2,3,4,5,6],fontsize = 15)
    ax[2].set_xlim(0,1)
    ax[2].set_xticks(np.round(np.arange(0,1.1,0.2),1))
    ax[2].set_xticklabels(np.round(np.arange(0,1.1,0.2),1), fontsize =15)

    os.makedirs("CFADs", exist_ok = True)
    plt.savefig("CFADs/lines_cfad"+id+".png", dpi = 300)
    plt.close()
