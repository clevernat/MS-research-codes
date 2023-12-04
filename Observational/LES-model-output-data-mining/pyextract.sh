#!/bin/bash

# Define the Python code as a multi-line string
python_code='''

print("Importing Libraries...")
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import xarray as xr

import cartopy.crs as ccrs
import cartopy.feature as cfeature

from datetime import datetime
from cartopy.util import add_cyclic_point
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import gc
import glob

from netCDF4 import Dataset

from matplotlib.colors import Normalize
import matplotlib.colors as mcolors
import matplotlib.patheffects as PathEffects
#======================
print("finished importing Libraries")

print("starting reading data")
path_to_variables = np.sort(glob.glob("/glade/campaign/ral/wsap/tjuliano/comble/coupled_meso_micro/RESTART_*/compressed/wrfout_cloud_d02_2020-03-*"))
path_to_lat_lon = np.sort(glob.glob("/glade/campaign/ral/wsap/tjuliano/comble/coupled_meso_micro/RESTART_*/compressed/wrfout_wind_d02_2020-03-*"))
print("done reading data")
#=================

print("Getting the height from the wind path")
# Assuming youve already opened the dataset and have the "HGT" variable in meters
wind_vars = xr.open_dataset(path_to_lat_lon[0])
height_meters = wind_vars["HGT"]  # "HGT" is in meters

# Convert height to kilometers
height_kilometers = height_meters / 1000
print("Done converting...")

# Initialize an empty list to store modified datasets
print("Loading function and performing computations")
datasets_list = []

for i, path in enumerate(path_to_variables[:20]):
    # List of variables you want to keep
    variables_to_keep = ["QRAIN", "QGRAUP", "QSNOW", "Times"]
    
    # Open the dataset
    ds = xr.open_dataset(path)
    
    # Remove all variables except the ones in variables_to_keep
    ds = ds.drop_vars([var for var in ds.data_vars if var not in variables_to_keep])
    
    # Add the height_below_1km variable to the dataset
    ds["HGT"] = height_kilometers  # Assuming the dimensions match
    
    # Filter the height variable to include values less than or equal to 1 km 
    ds = ds.where(ds["HGT"]<=1, drop=True)
    
    # Append the modified dataset to the list
    datasets_list.append(ds)
    
#==========
print("Done Loading function and performing computations")

print("Merging the three files")
# Combine the datasets in the list into a single dataset along the "Time" dimension
combined_dataset = xr.concat(datasets_list, dim="Time")
print("Done merging")

#======
print("Saving files to a new netcdf ...")
# Save the merged dataset to a new NetCDF file
output_file = "/glade/derecho/scratch/notengmerged_data.nc"  # Provide the desired output file name

merged_dataset.to_netcdf(output_file)
print("Done saving files")



'''

# Execute the Python code
echo "$python_code" | python -