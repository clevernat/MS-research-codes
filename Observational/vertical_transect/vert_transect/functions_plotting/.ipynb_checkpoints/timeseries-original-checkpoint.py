import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
from matplotlib.transforms import Bbox
from matplotlib.lines import Line2D

def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

def timeseries(vars, pn, bts, bte,  lwp_time, pwv_time, aeri, time, id, st):
    labels = [time[int((len(time)-1)*0.1*i)] for i in range(0,11)]
    bt = bte-bts
    if st[0:2] == "BI":
        ylim = int(np.round(np.nanmax(vars['kazr']['h'])/1000.0,0)*1000)
        yticks = range(0,1100,200)
        yticklabels = np.round(np.arange(0,1.1,0.2),1)
    if st[0:3] == "ANX":
        ylim = int(np.round(np.nanmax(vars['kazr']['h'])/1000.0,0)*1000)
        yticks = range(0,ylim+1,1000)
        yticklabels = range(0,int(ylim/1000+1),1)
    if st[0:3] == "NSA":
        ylim = int(np.round(np.nanmax(vars['kazr']['h'])/1000.0,0)*1000)
        yticks = range(0,2100,1000)
        yticklabels = range(0,3,1)

    lab = ["(a)","(b)","(c)","(d)","(e)","(f)", "(g)","(h)"]
    plottype = "cells"
    ticks = np.linspace(0,bt,num=11)

    ps = np.asarray(list(pn.values()))
    x = len(np.where(ps>=0)[0])

    #fig, ax = plt.subplots(x,1, figsize=(15,11/6*x))
    fig, ax = plt.subplots(x,1, figsize=(3,11/6*x))
    #fig, ax = plt.subplots(2,2, figsize=(3,11))
    if x == 1:
        ax = [ax]
        print(ax)

    if ps[0] >= 0:
        vmin = -30
        vmax = 25
        num = vmax - vmin + 1
        cmap = plt.get_cmap('nipy_spectral')
        bounds = np.linspace(vmin,vmax,num=num)
        colorlist  = [cmap(i) for i in range(cmap.N)]
        cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', colorlist, cmap.N)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

        circle = ax[ps[0]].pcolormesh(vars['kazr']['time']-bts, vars['kazr']['h'], np.transpose(vars['kazr']['ref']),cmap=cmap,vmin=vmin,vmax=vmax,zorder=5)
        #circlevpt = ax[ps[0]].pcolormesh(vars['vpt']['time']-bts, vars['vpt']['h'], np.transpose(vars['vpt']['ref']),cmap=cmap,vmin=vmin,vmax=vmax,zorder=4)
        ax[ps[0]].set_ylim(0,ylim)
        ax[ps[0]].set_xlim(0,bt)
        ax[ps[0]].set_yticks(yticks)
        #ax[ps[0]].set_yticklabels(yticklabels, fontsize=13)
        ax[ps[0]].set_yticklabels(yticklabels, fontsize=5) #new
        ax[ps[0]].set_xticks(ticks)
        #ax[ps[0]].set_ylabel("MSL [km]", fontsize=13)
        ax[ps[0]].set_ylabel("MSL [km]", fontsize=7.2) #new
        if ps[0] == max(ps):
            #ax[ps[0]].set_xticklabels(labels, fontsize=10.1)
            ax[ps[0]].set_xticklabels(labels, fontsize=2.3, fontweight='bold') #new
            ax[ps[0]].tick_params(labelbottom=True,labelleft=True, labelright=True,bottom=True, top=True, left=True, right=True)
        else:
            ax[ps[0]].tick_params(labelbottom=False,labelleft=True, labelright=True,bottom=True, top=True, left=True, right=True)
        ax[ps[0]].yaxis.set_label_position("right")
        ax[ps[0]].tick_params('both', length=6, width=1, which='major')
        #ax[ps[0]].set_title(lab[ps[0]],loc='left',fontsize=13,pad=3.0)
        
        cbar_ax = fig.add_axes([1.1, 0.11, 0.027, 0.76]) #lbwh #new


        #cbar=fig.colorbar(cs, cax=cbar_ax, orientation='vertical') #new
        cbar = fig.colorbar(circle, norm=norm, cax=cbar_ax, ax=ax[ps[0]],extend= "neither") #new
        #cbar = fig.colorbar(circle, norm=norm, ax=ax[ps[0]],extend= "neither")
        
        cbar.set_ticks([-30,-20,-10,0,10,20])
        cbar.set_ticklabels([-30,-20,-10,0,10,20])
        #cbar.ax.tick_params(labelsize=11)
        cbar.ax.tick_params(labelsize=7) # new
        #cbar.set_label("Reflectivity [dBZ]", fontsize=11.15)
        cbar.set_label("Reflectivity [dBZ]", fontsize=9.15) #new


    if ps[1] >=0:
        vmin = -3
        vmax = 1
        num = vmax*10 - vmin*10 + 1
        cmap = plt.get_cmap('bwr')
        bounds = np.linspace(vmin,vmax,num=num)
        colorlist  = [cmap(i) for i in range(cmap.N)]
        cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', colorlist, cmap.N)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

        circle = ax[ps[1]].pcolormesh(vars['kazr']['time']-bts, vars['kazr']['h'], np.transpose(vars['kazr']['vel']),cmap=cmap,vmin=vmin,vmax=vmax,zorder=5)
        #circlevpt = ax[ps[1]].pcolormesh(vars['vpt']['time']-bts, vars['vpt']['h'], np.transpose(vars['vpt']['vel']),cmap=cmap,vmin=vmin,vmax=vmax,zorder=4)
        ax[ps[1]].set_ylim(0,ylim)
        ax[ps[1]].set_xlim(0,bt)
        ax[ps[1]].set_ylabel("MSL [km]", fontsize=13)
        ax[ps[1]].set_yticks(yticks)
        ax[ps[1]].set_yticklabels(yticklabels, fontsize=13)
        ax[ps[1]].set_xticks(ticks)
        if ps[1] == max(ps):
            ax[ps[1]].set_xticklabels(labels, fontsize=10.1)
            ax[ps[1]].tick_params(labelbottom=True,labelleft=True, labelright=True,bottom=True, top=True, left=True, right=True)
        else:
            ax[ps[1]].tick_params(labelbottom=False,labelleft=True, labelright=True,bottom=True, top=True, left=True, right=True)
        ax[ps[1]].yaxis.set_label_position("right")
        ax[ps[1]].tick_params('both', length=6, width=1, which='major')
        ax[ps[1]].set_title(lab[ps[1]],loc='left',fontsize=13,pad=3.0)
        cbar = fig.colorbar(circle, norm=norm, ax=ax[ps[1]],extend="both")
        cbar.set_ticks(bounds[::10])
        cbar.set_ticklabels([-3,-2,-1,0,1])
        cbar.ax.tick_params(labelsize=11)
        cbar.set_label("Doppler\nVelocity [m s$^{-1}$]", fontsize=11.15)
    if ps[2] >= 0:
        vmin = 0
        vmax = 1
        num = vmax*10 - vmin*10 + 1
        cmap = plt.get_cmap('Blues')
        #cmap.set_under("k")
        bounds = np.linspace(vmin,vmax,num=num)
        colorlist  = [cmap(i) for i in range(cmap.N)]
        cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', colorlist, cmap.N)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
        circle = ax[ps[2]].pcolormesh(vars['kazr']['time']-bts, vars['kazr']['h'], np.transpose(vars['kazr']['sw']),cmap=cmap,vmin=vmin,vmax=vmax,zorder=5)
        #circlevpt = ax[ps[2]].pcolormesh(vars['vpt']['time']-bts, vars['vpt']['h'], np.transpose(vars['vpt']['sw']),cmap=cmap,vmin=vmin,vmax=vmax,zorder=4)
        ax[ps[2]].set_ylim(0,ylim)
        ax[ps[2]].set_xlim(0,bt)
        ax[ps[2]].set_ylabel("MSL [km]", fontsize=13)
        ax[ps[2]].set_yticks(yticks)
        ax[ps[2]].set_yticklabels(yticklabels, fontsize=13)
        ax[ps[2]].set_xticks(ticks)
        if ps[2] == max(ps):
            ax[ps[2]].set_xticklabels(labels, fontsize=10.1)
            ax[ps[2]].tick_params(labelbottom=True,labelleft=True, labelright=True,bottom=True, top=True, left=True, right=True)
        else:
            ax[ps[2]].tick_params(labelbottom=False,labelleft=True, labelright=True,bottom=True, top=True, left=True, right=True)
        ax[ps[2]].tick_params('both', length=6, width=1, which='major')
        ax[ps[2]].yaxis.set_label_position("right")
        ax[ps[2]].set_title(lab[ps[2]],loc='left',fontsize=13,pad=3.0)
        cbar = fig.colorbar(circle, norm=norm, ax=ax[ps[2]],extend= "max")
        cbar.set_ticks(np.round(bounds[::2],1))
        cbar.set_ticklabels(np.round(bounds[::2],1))
        cbar.ax.tick_params(labelsize=11)
        cbar.set_label("Spectral\nWidth [m s$^{-1}$]", fontsize=11.75)

    if ps[3] >= 0:
        # vmin = 0
        # vmax = 0.5
        # num = vmax*10 - vmin*10 + 1
        # cmap = plt.get_cmap('jet')
        # cmap.set_under("k")
        # bounds = np.linspace(vmin,vmax,num=num)
        # colorlist  = [cmap(i) for i in range(cmap.N)]
        # cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', colorlist, cmap.N)
        # norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

        color = plt.cm.get_cmap('viridis')
        colorlist  = [color(k) for k in range(color.N)]
        #color_tup = [(256/256, 256/256, 256/256, 1.0), (0/256, 0/256, 200/256, 1.0), (152/256, 245/256, 255/256, 1.0), (255/256, 0/256, 0/256, 1.0), (171/256, 130/256, 255/256, 1.0), (0.0,1.0,0.0,1.0),(255/256, 128/256, 0, 1.0), (0,0,0,1.0), (127/256, 127/256, 127/256, 1.0)]
        color_tup = [(256/256, 256/256, 256/256, 1.0), (0/256, 0/256, 200/256, 1.0), (152/256, 245/256, 255/256, 1.0), (255/256, 0/256, 0/256, 1.0), (0,0,0,1.0)]
        # for k1 in range(0,9):
        #     for k2 in range(int(len(colorlist)/9*k1), int(len(colorlist)/9*k1+len(colorlist)/9)):
        #         colorlist[k2] = color_tup[k1]
        for k1 in range(0,5):
             for k2 in range(int(len(colorlist)/5*k1), int(len(colorlist)/5*k1+len(colorlist)/5)):
                 colorlist[k2] = color_tup[k1]
        color = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', colorlist, color.N)
        #bounds = [0,1,2,3,4,5,6,7,8,9]
        bounds = [0,1,2,3,4,8]
        norm = mpl.colors.BoundaryNorm(bounds, color.N)

        circle = ax[ps[3]].pcolormesh(vars['cld']['time']-bts, vars['cld']['h'], np.transpose(vars['cld']['phase']),cmap=color,norm=norm,vmin=0,vmax=9)

        ax[ps[3]].set_ylim(0,ylim)
        ax[ps[3]].set_ylabel("MSL [km]", fontsize=13)
        ax[ps[3]].set_yticks(yticks)
        ax[ps[3]].set_yticklabels(yticklabels, fontsize=13)
        ax[ps[3]].set_xticks(ticks)
        if ps[3] == max(ps):
            ax[ps[3]].set_xticklabels(labels, fontsize=10.1)
            ax[ps[3]].tick_params(labelbottom=True,labelleft=True, labelright=True,bottom=True, top=True, left=True, right=True)
        else:
            ax[ps[3]].tick_params(labelbottom=False,labelleft=True, labelright=True,bottom=True, top=True, left=True, right=True)
        ax[ps[3]].yaxis.set_label_position("right")
        ax[ps[3]].tick_params('both', length=6, width=1, which='major')
        ax[ps[3]].set_title(lab[ps[3]],loc='left',fontsize=13,pad=3.0)
        cbar = fig.colorbar(circle, norm=norm, ax=ax[ps[3]],extend= "neither")
        # cbar.set_ticks(np.round(bounds[::2],1))
        # cbar.set_ticklabels(np.round(bounds[::2],1))
        #cbar.set_ticks([0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5])
        #cbar.set_ticklabels(["Clear", "Liq", "Ice", "Mixed", "Driz", "L+D","Rain","Snow","Unkn"])
        cbar.set_ticks([0.5,1.5,2.5,3.5,6])
        cbar.set_ticklabels(["Clear", "Liq", "Ice", "Mixed","Snow"])
        cbar.ax.tick_params(labelsize=11)
        cbar.set_label("Cloud Phase", fontsize=12)

    if ps[4] >= 0:
        axpwv = ax[ps[4]].twinx()
        axpwv.spines['left'].set_color('blue')
        axpwv.spines['right'].set_color('red')
        #axiwp = ax[ps[4]].twinx()
        #axiwp.spines["right"].set_position(("axes", 1.08))
        #make_patch_spines_invisible(axiwp)
        #axiwp.spines["right"].set_visible(True)
        ax[ps[4]].scatter(lwp_time-bts, vars['mwr']['lwp'], c = 'blue', s=1.225, zorder = 2)
        nn = axpwv.scatter(pwv_time[:740]-bts, vars['mwr']['pwv'][:740]*10, c = 'red', s=1.225, zorder = 2)
        #axiwp.scatter(vars['cloud']['time'][::8]-bts,vars['iwp']['iwp'][::8], c = 'k', s=1.225, alpha = 0.8, zorder = 2)
        #axiwp.set_ylabel("IWP [kg m$^{-2}$]", c = "k", fontsize=10.5)
        axpwv.set_ylabel("PWV [kg m$^{-2}$]", c = "red", fontsize=13)
        ax[ps[4]].set_ylabel("LWP [kg m$^{-2}$]", c = "blue", fontsize=13)
        ax[ps[4]].tick_params(axis='y', colors='blue')
        axpwv.tick_params(axis='y', colors='red')
        ax[ps[4]].set_xticks(ticks)
        ax[ps[4]].set_title(lab[ps[4]],loc='left',fontsize=13,pad=3.0)
        if len(vars['mwr']['lwp'])>0:
            lwp_ylim = int(np.round((np.nanmax(vars['mwr']['lwp'])+100)/100.0,0)*100)
            lwp_yticks = int(np.round(lwp_ylim/100.0,0)*20)
            #iwp_ylim = int(np.round((np.nanmax(vars['iwp']['iwp'])+100)/100.0,0)*100)
            #iwp_yticks = int(np.round(iwp_ylim/100.0,0)*20)
            for i in range(1,6):
                pwv_ylim = int(np.round((np.nanmax(vars['mwr']['pwv'][:740]*10)+i),0))
                if pwv_ylim % 5 == 0:
                    break
            pwv_yticks = int(np.round(pwv_ylim,0)*0.2)
            if np.isnan(np.nanmax(vars['mwr']['pwv'][:740]*10)) == True :
                pwv_ylim = 100
            axpwv.set_ylim(0,pwv_ylim)
            axpwv.set_yticks(range(0,pwv_ylim+1,pwv_yticks))
            axpwv.set_yticklabels(np.round(np.arange(0,pwv_ylim+1,pwv_yticks),0), fontsize=13)
            ax[ps[4]].set_ylim(0,lwp_ylim)
            #ax[ps[4]].set_yticks(range(0,lwp_ylim+1,lwp_yticks))
            #ax[ps[4]].set_yticklabels(np.round(np.arange(0,lwp_ylim/1000+1,lwp_yticks/1000),2), fontsize=13)
            #axiwp.set_ylim(0,iwp_ylim)
            #axiwp.set_yticks(range(0,iwp_ylim+1,iwp_yticks))
            #axiwp.set_yticklabels(np.round(np.arange(0,iwp_ylim/1000+1,iwp_yticks/1000),2), fontsize=13)
        ax[ps[4]].set_xlim(0,bt)
        if ps[4] == max(ps):
            ax[ps[4]].set_xticklabels(labels, fontsize=10.1)
            ax[ps[4]].tick_params(labelbottom=True,labelleft=True, labelright=False,bottom=True, top=True, left=True, right=False)
        else:
            ax[ps[4]].tick_params(labelbottom=False,labelleft=True, labelright=False,bottom=True, top=True, left=True, right=False)
        ax[ps[4]].tick_params('both', length=6, width=1, which='major')
        ax[ps[4]].grid(axis="y", ls = '--', c = "grey",alpha = 0.5)
        cbar_ps4 = fig.colorbar(nn, ax=ax[ps[4]],extend= "neither")
        #cbar_ps4 = fig.colorbar(circle, ax=ax[ps[4]],extend= "neither")
        cbar_ps4.set_ticks(np.round(bounds[::2],1))
        cbar_ps4.set_ticklabels(np.round(bounds[::2],1))
        cbar_ps4.ax.tick_params(labelsize=10)
        cbar_ps4.set_label("LDR", fontsize=8)
        cbar_ps4.remove()

        ### Precipitation ###
        ax7t = ax[ps[4]].twinx()
        ax7t.spines["right"].set_position(("axes", 1.08))
        make_patch_spines_invisible(ax7t)
        ax7t.spines["right"].set_visible(True)
        #circle3 = ax7t.plot(vars['sur']['time']-bts,np.nancumsum(vars['sur']['prec']/60), c = 'k',linewidth=1.5)
        circle3 = ax7t.scatter(vars['sur']['time']-bts,vars['sur']['prec'], c = 'black', marker='+', s=1.225, zorder = 2)
        for i in range(1,3):
            prec_ylim = int(np.round((max(vars['sur']['prec'])+i-0.5),0))
            if prec_ylim % 2 == 0:
                break
        ax7t.set_ylim(0,prec_ylim)
        ax7t.set_yticks(np.round(np.linspace(0,prec_ylim,5),1))
        ax7t.set_yticklabels(np.round(np.linspace(0,prec_ylim,5),1), fontsize=13)
        ax7t.set_ylabel("Precipitation\nRate [mm hr$^{-1}$]", c = "k", fontsize=12)

    if ps[5] >= 0:
        db = int(len(vars['sur']['wspd'])/30.0)
        if db == 0:
            db = 1
        ax6t = ax[ps[5]].twinx()
        ax6t.spines['left'].set_color('red')
        ax6t.spines['right'].set_color('blue')

        ### Precipitation ###
        # ax7t = ax[ps[5]].twinx()
        # ax7t.spines["right"].set_position(("axes", 1.08))
        # make_patch_spines_invisible(ax7t)
        # ax7t.spines["right"].set_visible(True)
        # circle3 = ax7t.plot(vars['sur']['time']-bts,np.nancumsum(vars['sur']['prec']/60), c = 'k',linewidth=1.5)
        # for i in range(1,3):
        #     prec_ylim = int(np.round((max(np.nancumsum(vars['sur']['prec']/60))+i-0.5),0))
        #     if prec_ylim % 2 == 0:
        #         break
        # ax7t.set_ylim(0,prec_ylim)
        # ax7t.set_yticks(np.round(np.linspace(0,prec_ylim,5),1))
        # ax7t.set_yticklabels(np.round(np.linspace(0,prec_ylim,5),1), fontsize=13)
        # ax7t.set_ylabel("Precipitation [mm]", c = "k", fontsize=10.5)

        ax[ps[5]].set_title(lab[ps[5]],loc='left',fontsize=13,pad=3.0)
        maxax6, minax6 = int(np.nanmax(vars['sur']['t'])+1.2), int(np.nanmin(vars['sur']['t'])-1.0)
        maxax6t, minax6t = int(np.nanmax(vars['sur']['pres'])+1.2), int(np.nanmin(vars['sur']['pres'])-0.2)

        circle1 = ax[ps[5]].plot(vars['sur']['time']-bts,vars['sur']['t'], c = 'red',linewidth=1.5)
        circle2 = ax6t.plot(vars['sur']['time']-bts,vars['sur']['pres'], c = 'blue',linewidth=1.5)

        ax[ps[5]].barbs(vars['sur']['time'][int(db/2)::db]-bts,np.full(len(vars['sur']['t']),minax6+(maxax6-minax6)/2)[int(db/2)::db],-1.944*vars['sur']['wspd'][int(db/2)::db]*np.sin(vars['sur']['dirc'][int(db/2)::db]*np.pi/180.0),-1.944*vars['sur']['wspd'][int(db/2)::db]*np.cos(vars['sur']['dirc'][int(db/2)::db]*np.pi/180.0),pivot='middle',color='black', length = 7)
        ax[ps[5]].set_ylabel("T [°C]", c = "red", fontsize=13)
        ax[ps[5]].set_xlim(0,bt)
        ax[ps[5]].set_ylim(minax6,maxax6)
        ax6t.set_ylim(minax6t,maxax6t)
        ax[ps[5]].set_xticks(ticks)
        if ps[5] == max(ps):
            #ax[ps[5]].set_xticklabels(labels, fontsize=10.1)
            ax[ps[5]].set_xticklabels(labels, fontsize=2.1)
            ax[ps[5]].tick_params(labelbottom=True,labelleft=True, labelright=False,bottom=True, top=True, left=True, right=False)
        else:
            ax[ps[5]].tick_params(labelbottom=False,labelleft=True, labelright=False,bottom=True, top=True, left=True, right=False)
        ax[ps[5]].tick_params(axis='y', colors='red', labelsize=13)
        ax6t.tick_params(axis='y', colors='blue', labelsize=13)
        ax6t.set_ylabel("P [hPa]", c = "blue", fontsize=13)
        ax6t.tick_params('both', length=6, width=1, which='major')
        ax[ps[5]].grid(axis="y", ls = '--', c = "red",alpha = 0.45)
        ax[ps[5]].tick_params('both', length=6, width=1, which='major')
        ax6t.grid(axis="y", ls = '--', c = "blue",alpha = 0.45)
        cbar = fig.colorbar(nn, ax=ax[ps[5]])
        cbar.remove()

    if ps[6] >= 0:
        vmin = -2
        vmax = 2
        num = vmax*10 - vmin*10 + 1
        cmap = plt.get_cmap('bwr')
        cmap.set_under("k")
        bounds = np.linspace(vmin,vmax,num=num)
        colorlist  = [cmap(i) for i in range(cmap.N)]
        cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', colorlist, cmap.N)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
        circle = ax[ps[6]].pcolormesh(aeri[2], vars['aeri']['h'], np.transpose(aeri[0]),cmap=cmap,norm=norm,vmin=vmin,vmax=vmax)
        cntr2 = ax[ps[6]].contour(aeri[2], vars['aeri']['h'], np.transpose(aeri[1]),np.round(np.arange(0,150,10)), colors='k', linewidths=0.8, linestyles ="dashed", zorder=6)
        plt.clabel(cntr2,inline=1, inline_spacing=8, fontsize=7, fmt='%i')
        ax[ps[6]].set_ylim(0,3000)
        ax[ps[6]].set_xlim(0,bt)
        ax[ps[6]].set_ylabel("MSL [km]", fontsize=13)
        ax[ps[6]].set_yticks([0,1000,2000,3000])
        ax[ps[6]].set_yticklabels([0,1,2,3], fontsize=13)
        ax[ps[6]].set_xticks(ticks)
        if ps[6] == max(ps):
            ax[ps[6]].set_xticklabels(labels, fontsize=10.1)
            ax[ps[6]].tick_params(labelbottom=True,labelleft=True, labelright=True,bottom=True, top=True, left=True, right=True)
        else:
            ax[ps[6]].tick_params(labelbottom=False,labelleft=True, labelright=True,bottom=True, top=True, left=True, right=True)
        ax[ps[6]].tick_params('both', length=6, width=1, which='major')
        ax[ps[6]].yaxis.set_label_position("right")
        ax[ps[6]].set_title(lab[ps[6]],loc='left',fontsize=13,pad=3.0)
        cbar = fig.colorbar(circle, norm=norm, ax=ax[ps[6]],extend= "both")
        cbar.set_ticks([-2,-1,0,1,2])
        cbar.set_ticklabels([-2,-1,0,1,2])
        cbar.ax.tick_params(labelsize=11)
        cbar.set_label("Temperature \n Perturbation [°C]", fontsize=12)

    if ps[8] >= 0:
        vmin = -5
        vmax = 1
        num = vmax - vmin
        cmap = plt.get_cmap('jet')
        bounds = np.linspace(-5,1,num=31)
        cmap = plt.get_cmap('jet')
        colorlist  = [cmap(i) for i in range(cmap.N)]
        cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', colorlist, cmap.N)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
        invalid_iwc = np.greater(vars['micro']['iwc'],1)
        vars['micro']['iwc'][invalid_iwc]=np.nan
        circle = ax[ps[8]].pcolormesh(vars['micro']['time']-bts, vars['micro']['h'], np.transpose(np.log10(vars['micro']['iwc'])),cmap=cmap,norm=norm,vmin=vmin,vmax=vmax)
        ax[ps[8]].set_ylim(0,ylim)
        ax[ps[8]].set_xlim(0,bt)
        ax[ps[8]].set_yticks(yticks)
        ax[ps[8]].set_yticklabels(yticklabels, fontsize=13)
        ax[ps[8]].set_xticks(ticks)
        ax[ps[8]].set_ylabel("MSL [km]", fontsize=13)
        if ps[8] == max(ps):
            ax[ps[8]].set_xticklabels(labels, fontsize=10.1)
            ax[ps[8]].tick_params(labelbottom=True,labelleft=True, labelright=True,bottom=True, top=True, left=True, right=True)
        else:
            ax[ps[8]].tick_params(labelbottom=False,labelleft=True, labelright=True,bottom=True, top=True, left=True, right=True)
        ax[ps[8]].yaxis.set_label_position("right")
        ax[ps[8]].tick_params('both', length=6, width=1, which='major')
        ax[ps[8]].set_title(lab[ps[8]],loc='left',fontsize=13,pad=3.0)
        cbar = fig.colorbar(circle, norm=norm, ax=ax[ps[8]],extend= "neither")
        cbar.set_ticks(bounds[::10])
        cbar.set_ticklabels(bounds[::10])
        cbar.ax.tick_params(labelsize=10)
        cbar.set_label("log$_{10}$(IWC) [g m$^{-3}$]", fontsize=10.5)

    # ax[max(ps)+1].set_xlabel("Time [UTC]", fontsize=13)
    # ax[max(ps)+1].tick_params(labelbottom=True,labelleft=False, labelright=False,bottom=True, top=True, left=False, right=False)
    # ax[max(ps)+1].set_xlim(0,bt)
    # ax[max(ps)+1].set_ylim(0,1)
    # ax[max(ps)+1].set_xticks(ticks)
    # ax[max(ps)+1].set_xticklabels(labels, fontsize=10.1)
    # ax[max(ps)+1].tick_params('both', length=6, width=1, which='major')
    # ax[max(ps)+1].scatter(30420,0.98,c="k", marker="v",s=125)
    # ax[max(ps)+1].axvline(6300,0,1,c="blue", linewidth=1.25, zorder=6, linestyle="dashed")
    # ax[max(ps)+1].axvline(19800,0,1,c="blue", linewidth=1.25, zorder=6, linestyle="dashed")
    # ax[max(ps)+1].axvline(28800,0,1,c="blue", linewidth=1.25, zorder=6, linestyle="dashed")
    # ax[max(ps)+1].axvline(16200,0,1,c="k", linewidth=1.75, zorder=6, linestyle="solid")
    # ax[max(ps)+1].axvline(24000,0,1,c="k", linewidth=1.75, zorder=6, linestyle="solid")
    # ax[max(ps)+1].scatter(15540,0.98,c="k", marker="v",s=125)
    # ax[max(ps)+1].axvline(16200,0,ylim,c="blue", linewidth=1.25, zorder=6, linestyle="dashed")
    # ax[max(ps)+1].axvline(19800,0,ylim,c="blue", linewidth=1.25, zorder=6, linestyle="dashed")
    # ax[max(ps)+1].axvline(25200,0,ylim,c="blue", linewidth=1.25, zorder=6, linestyle="dashed")
    # ax[max(ps)+1].axvline(10800,0,1,c="k", linewidth=1.75, zorder=6, linestyle="solid")
    # ax[max(ps)+1].axvline(18300,0,1,c="k", linewidth=1.75, zorder=6, linestyle="solid")
    # ax[max(ps)+1].axvline(21600,0,1,c="k", linewidth=1.75, zorder=6, linestyle="solid")
    # ax[max(ps)+1].axvline(28800,0,1,c="k", linewidth=1.75, zorder=6, linestyle="solid")
    # ax[max(ps)+1].scatter(10620,0.98,c="k", marker="v",s=125)
    # ax[max(ps)+1].axvline(12600,0,ylim,c ="blue", linewidth=1.25, zorder=6, linestyle="dashed")
    # ax[max(ps)+1].axvline(18900,0,ylim,c ="blue", linewidth=1.25, zorder=6, linestyle="dashed")
    # ax[max(ps)+1].axvline(25200,0,ylim,c ="blue", linewidth=1.25, zorder=6, linestyle="dashed")
    # ax[max(ps)+1].axvline(5400,0,1,c="k", linewidth=1.75, zorder=6, linestyle="solid")
    # ax[max(ps)+1].axvline(9600,0,1,c="k", linewidth=1.75, zorder=6, linestyle="solid")
    # ax[max(ps)+1].axvline(14580,0,1,c="k", linewidth=1.75, zorder=6, linestyle="solid")
    # ax[max(ps)+1].axvline(21600,0,1,c="k", linewidth=1.75, zorder=6, linestyle="solid")
    # ax[max(ps)+1].axvline(21600,0,1,c="k", linewidth=1.75, zorder=6, linestyle="solid")
    # ax[max(ps)+1].axvline(30600,0,1,c="k", linewidth=1.75, zorder=6, linestyle="solid")

    # ax[max(ps)+1].set_title(lab[max(ps)+1],loc='left',fontsize=13,pad=3.0)
    # cbar = fig.colorbar(circle, norm=norm, ax=ax[max(ps)+1])
    # cbar.remove()
    if ps[0] >= 0:
        ### Vertical Temperature ###
        bin = np.linspace(0, len(vars['is']['t']),len(vars['is']['t']))
        y = vars['is']['h']
        cntr = ax[ps[0]].contour(vars['is']['time'][::2]-bts,y[::2],np.transpose(vars['is']['t'])[::2,::2],np.round(np.arange(-100,320,10)), colors='k', linewidths=0.8, linestyles ="dashed", zorder=6)
        plt.clabel(cntr,inline=1, inline_spacing=8, fontsize=8, fmt='%i',manual=False)

    os.makedirs("t_series/t_series_"+st, exist_ok = True)
    plt.savefig("t_series/t_series_"+st+"/t_series_"+st+"_"+plottype+"_"+id+".png", dpi = 330, bbox_inches='tight')
    plt.close()
