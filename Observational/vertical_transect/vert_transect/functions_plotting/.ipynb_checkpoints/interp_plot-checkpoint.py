import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
from matplotlib.transforms import Bbox
from matplotlib.lines import Line2D

def interp_plot(interp,wind,wind2,height,time,id,st):
    ticks = np.linspace(0,len(interp[0]),num=11, dtype = "int")
    labels = [time[int((len(time)-1)*0.1*i)] for i in range(0,11)]
    yticks = range(0,7100,1000)
    yticklabels = range(0,8,1)

    x = len(interp)
    fig, ax = plt.subplots(x,1,figsize = (10.9,10/6*x))

    vmin = int(min(np.ravel(interp[0]))-1.0)
    vmax = int(max(np.ravel(interp[0]))+1.0)
    num = vmax - vmin + 1
    cmap = plt.get_cmap('jet')
    cmap.set_under("k")
    bounds = np.linspace(vmin,vmax,num=num)
    colorlist  = [cmap(i) for i in range(cmap.N)]
    cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', colorlist, cmap.N)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    bin = np.linspace(0, len(interp[0]),len(interp[0])+1)

    circle = ax[0].pcolormesh(bin, height, np.transpose(interp[0]),cmap=cmap,norm=norm,vmin=vmin,vmax=vmax)
    ax[0].set_ylim(0,7000)
    ax[0].set_ylabel("MSL [km]")
    ax[0].set_yticks(yticks)
    ax[0].set_yticklabels(yticklabels)
    ax[0].set_xticks(ticks)
    ax[0].tick_params(labelbottom=False,labelleft=True, labelright=True,bottom=True, top=True, left=True, right=True)
    ax[0].yaxis.set_label_position("right")
    cbar = fig.colorbar(circle, norm=norm, ax=ax[0],extend= "neither")
    cbar.set_ticks(np.round(bounds[::5],1))
    cbar.set_ticklabels(np.round(bounds[::5],1))
    cbar.ax.tick_params(labelsize=10)
    cbar.set_label("$\\theta_{e}$ [K]", fontsize=8)

    v_list = [abs(int(min(np.ravel(interp[1]))-1.0)),abs(int(max(np.ravel(interp[1]))+1.0))]
    vmin = -1*max(v_list)
    vmax = max(v_list)
    num = (vmax - vmin)*10 + 1
    cmap = plt.get_cmap('bwr')
    cmap.set_under("k")
    bounds = np.linspace(vmin,vmax,num=num)
    colorlist  = [cmap(i) for i in range(cmap.N)]
    cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', colorlist, cmap.N)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    circle = ax[1].pcolormesh(bin, height, np.transpose(interp[1]),cmap=cmap,norm=norm,vmin=vmin,vmax=vmax)
    ax[1].set_ylim(0,7000)
    ax[1].set_ylabel("MSL [km]")
    ax[1].set_yticks(yticks)
    ax[1].set_yticklabels(yticklabels)
    ax[1].set_xticks(ticks)
    ax[1].tick_params(labelbottom=False,labelleft=True, labelright=True,bottom=True, top=True, left=True, right=True)
    ax[1].yaxis.set_label_position("right")
    cbar = fig.colorbar(circle, norm=norm, ax=ax[1],extend= "neither")
    cbar.set_ticks(np.round(bounds[::10],1))
    cbar.set_ticklabels(np.round(bounds[::10],1))
    cbar.ax.tick_params(labelsize=10)
    cbar.set_label("Temperature \n Anomaly [K]", fontsize=8)

    v_list = [abs(int(min(np.ravel(interp[2]))-1.0)),abs(int(max(np.ravel(interp[2]))+1.0))]
    vmin = -1*max(v_list)
    vmax = max(v_list)
    num = (vmax - vmin) + 1
    cmap = plt.get_cmap('bwr')
    cmap.set_under("k")
    bounds = np.linspace(vmin,vmax,num=num)
    colorlist  = [cmap(i) for i in range(cmap.N)]
    cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', colorlist, cmap.N)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    circle = ax[2].pcolormesh(bin, height, np.transpose(interp[2]),cmap=cmap,norm=norm,vmin=vmin,vmax=vmax)
    ax[2].set_ylim(0,7000)
    ax[2].set_xlim(0,len(interp[2]))
    ax[2].set_ylabel("MSL [km]")
    ax[2].set_yticks(yticks)
    ax[2].set_yticklabels(yticklabels)
    ax[2].set_xticks(ticks)
    ax[2].set_xticklabels(labels)
    ax[2].tick_params(labelbottom=True,labelleft=True, labelright=True,bottom=True, top=True, left=True, right=True)
    ax[2].yaxis.set_label_position("right")
    if len(wind[0])>0:
        for i in range(0,len(wind[2]))[::2]:
            ax[2].barbs(wind[3],np.full(len(wind[3]),wind[2][i]),1.944*wind[0][:,i],1.944*wind[1][:,i],pivot='middle',color='black', length = 4)
        for i in range(0,len(wind2[2]))[::6]:
            ax[2].barbs(wind2[3],np.full(len(wind2[3]),wind2[2][i]),1.944*wind2[0][:,i],1.944*wind2[1][:,i],pivot='middle',color='black', length = 4)
    cbar = fig.colorbar(circle, norm=norm, ax=ax[2],extend= "neither")
    cbar.set_ticks(np.round(bounds[::10],1))
    cbar.set_ticklabels(np.round(bounds[::10],1))
    cbar.ax.tick_params(labelsize=10)
    cbar.set_label("Meridional Wind [m s$^{-1}$]", fontsize=8)

    os.makedirs("t_series/t_series_interp", exist_ok = True)
    plt.savefig("t_series/t_series_interp/t_series_interp_"+st+"_"+id+".png", dpi = 300)
    plt.close()
