import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
from matplotlib.transforms import Bbox
from matplotlib.lines import Line2D

# cfad(radar data, bins (x-axis), bins (y-axis) (radar range gates), len(bins(y-axis)), type of data (ref, vel, sw), id for file name,
# bins x-axis e.g.: (np.linspace(min,max,num=61)) from -30 dBZ to +30 dBZ with 61 bins
def cfad(data ,bin, height, h, type, id):
    invalid = np.logical_or(data>max(bin), data<min(bin))
    data[invalid] = np.nan
    ### Initalize arrays ###
    frequency = np.zeros((len(bin)-1,h),dtype=int) # data shown
    ref_ave = np.zeros(h)
    ref_90 = np.zeros(h)
    ref_10 = np.zeros(h)
    data_trans = np.transpose(data)

    ### Loop over every y bin ###
    for j in range(0, h):
        ### mean and 10th and 90th percentile for reflectivity ###
        if type == "ref":
            ref_ave[j] = np.nanmean(10**(data_trans[j]/10.0))
            ref_10[j] = np.nanpercentile(data_trans[j],10)
            ref_90[j] = np.nanpercentile(data_trans[j],90)
            ref_ave[j] = 10.0*np.log10(ref_ave[j])
        ### mean and 10th and 90th percentile for doppler velocity and spectrum width ###
        if type == "vel" or type == "sw":
            ref_ave[j] = np.nanmean(data_trans[j])
            ref_10[j] = np.nanpercentile(data_trans[j],10)
            ref_90[j] = np.nanpercentile(data_trans[j],90)
        ### loop over every x bin and assign data to each CFAD pixel ###
        for k in range(0,len(bin)-1):
            freq = np.where(np.logical_and(data_trans[j] >= bin[k], data_trans[j] < bin[k+1]))[0]
            frequency[k,j] = len(freq)

    ### Normalized frequency ###
    sum_frequency = np.sum(frequency)
    frequency = frequency/sum_frequency

    ### Minimum and Maximum frequecny for colorbar ###
    min_scale = np.nanmin(list(map(min, frequency)))
    max_scale = np.nanmax(list(map(max, frequency)))

    ### Colorbar parameters ###
    cmap = plt.get_cmap('viridis')
    bounds = np.arange(0,max_scale+0.00001,0.00001)
    colorlist  = [cmap(i) for i in range(cmap.N)]
    cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', colorlist, cmap.N)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    ### Plot CFAD ###
    fig, ax = plt.subplots(1,1, figsize=(10,6))
    cf = ax.pcolormesh(bin, height, np.transpose(frequency),cmap=cmap,norm=norm,vmin=min_scale,vmax=max_scale)

    ### Line plots of mean, 10th and 90th percentile ####
    hlim = 205
    ax.plot(ref_ave[:hlim], height[:hlim], c = 'red',linewidth=2.3)
    ax.plot(ref_10[:hlim], height[:hlim], ls = "--", c = 'red', linewidth=2.3)
    ax.plot(ref_90[:hlim], height[:hlim], ls = "dotted", c = 'red',linewidth=2.3)

    ### Differences between radar variables ###
    if type == "ref":
        ax.set_xlabel("Reflectivity [dBZ]",fontsize = 15)
        name = "_reflectivity"

        xticks = np.arange(-30,31,10)
    if type == "vel":
        ax.set_xlabel("Doppler Velocity [m s$^{-1}$]",fontsize = 15)
        name = "_velocity"
        cbx = 10
        ax.vlines(0,0,10000,linestyles="--",colors="white") # vertical lineat 0 m/s
        xticks = np.round(np.arange(-2.5,1.6,0.5),1)
    if type == "sw":
        ax.set_xlabel("Spectrum Width [m s$^{-1}$]",fontsize = 15)
        name = "_spectrum_width"
        cbx = 10
        xticks = np.round(np.arange(0,1.1,0.1),1)

    ### Labeling ###
    ax.set_ylabel("MSL [km]",fontsize = 15)
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.set_yticks(np.arange(0,10000,1000))
    ax.set_yticklabels(np.arange(0,10,1),fontsize = 15)
    ax.set_ylim(0,6000)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticks, fontsize =15)

    ### add colorbar to figure ###
    cbx = int(len(bounds)/10) # skip how many ticks of colorbar
    fig.subplots_adjust(right=0.8)
    cax = fig.add_axes([0.825,0.15,0.0225,0.7])
    cbar = fig.colorbar(cf, cax=cax)
    cbar.set_ticks(bounds[::cbx])
    cbar.set_ticklabels(np.around(bounds[::cbx]*10000, decimals=1))
    cbar.ax.tick_params(labelsize=14)
    cbar.set_label("Normalized Frequency [10$^{-4}$]", fontsize = 15)

    ### Save CFAD ###
    os.makedirs("CFAD", exist_ok = True)
    plt.savefig("CFAD/cfad"+name+id+".png", dpi = 300, bbox_inches='tight')
    plt.close()
