### Import modules required for this script ###
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import cartopy
import cartopy.feature as cfeature
import netCDF4
from datetime import datetime

### file path of radar file, needs to be changed to the file path you are using
path = "/glade/work/clackner/radar_data_for_nathaniel/"
### name of the file ###
file = "yrwms-nordic.mos.pcappi-0-dbz.noclass-clfilter-novpr-clcorr-block.laea-yrwms-1000.20200313.nc"

### Import file, to see more information use ncdump -h command on file in its directory ###
f = netCDF4.Dataset(path+file) # use netCDF4 library to read file
x = np.array(f.variables['lon'][:]) # import longitude, i.e. x axis
y = np.array(f.variables['lat'][:]) # import latitude, i.e. y axis
time = np.array(f.variables['time'][:]) # import time, here in base time
ref = np.array(f.variables['equivalent_reflectivity_factor'][:]) # import Reflectivity
f.close() # close opened file



### for all entries in time convert time value (base time) to desired format ###
time_title = [] # initalize list
for t in time:
    dt_object = datetime.utcfromtimestamp(t) # create datetime object
    time_title.append(dt_object.strftime('%m/%d/%Y %H:%M')) # append datetime object to list in desired information

### specify range in which to display radar ###
invalid = np.logical_or(ref > 40, ref < -20) # find values smaller than -20 or larger than 40
ref[invalid] = np.nan # set those values to nanmax

### specify color bar ###
vmin = -20 # min value of colorbar
vmax = 40 # max value of colorbar
color = plt.cm.get_cmap('jet') # colormap of colorbar can be changed to your liking
colorlist  = [color(i) for i in range(color.N)]
color = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', colorlist, color.N)
bounds = range(vmin,vmax+1) # array how divide colorbar
norm = mpl.colors.BoundaryNorm(bounds, color.N)

### center points and extent of map ###
latA = 69.141281 #latitude of COMBLE site
lonA = 15.684166-1 #longitude of COMBLE site -1
xm, ym = 6,2.75 # extent from center point in lon and lat
i = 0 # chose time step to plot

### Create figure ###
fig, ax = plt.subplots(1,1,figsize=(6,6),subplot_kw={'projection': ccrs.Orthographic(lonA,latA)}) # initalize figure with Orthographic projection
extent = [lonA-xm, lonA+xm, latA-ym, latA+ym] # define extent map
ax.set_extent(extent) # set extent of map
ax.coastlines(resolution='10m') # plot coastlines with high resolution: 10m

colormesh = ax.pcolormesh(x, y, ref[i], cmap=color, vmin=vmin, vmax=vmax,transform=ccrs.PlateCarree()) # plot Reflectivity pixels
ax.plot(lonA+1,latA, color='red', marker='*', markersize = 7.5,transform=ccrs.PlateCarree()) # plot red star at location of Andenes

fig.subplots_adjust(right=0.82) # create space for colorbar
cax = fig.add_axes([0.84,0.15,0.025,0.7]) # create new axes on which colorbar is plotted

### Color bar ###
cbar = fig.colorbar(colormesh,cax=cax)
cbar.set_ticks(bounds[::10])
cbar.set_ticklabels(bounds[::10])
cbar.ax.tick_params(labelsize=14)
cbar.set_label("Reflectivity Factor [dBZ]", fontsize = 15)

ax.set_title("Nordic Radar Mosaic: " +time_title[i]+" UTC") # title of figure

### parameter to name radar image based on time step i ###
stri = str(i)
if i < 100:
    stri = "0"+stri
    if i < 10:
        stri = "0"+stri

### Save and close figure ###
plt.savefig(path+"2020073_"+stri,dpi=300)
plt.close()
