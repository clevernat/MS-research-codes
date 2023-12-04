import time  # Import the time module for time-related functions
import glob  # Import the glob module for file path manipulation
import pandas as pd  # Import pandas for data manipulation
import xarray as xr  # Import xarray for working with NetCDF files
import gc  # Import the garbage collection module
import os  # Import the os module for file operations



# Get the current time to measure overall script execution time
start_time = time.time()

# Define the path to the NetCDF files
path = sorted(glob.glob("/glade/u/home/noteng/noteng/LES-output/wrfout_cloud_d02_2020-03-*"))

# Sample time_stamps
time_stamps = pd.date_range(start="2020-03-13 00:00:00", periods=len(path), freq="5T")  # Adjusted for the length of 'path'

# Initialize an empty list to store processed datasets
dataset = []

# Loop through all files in the list
for i, file_path in enumerate(path):
    start_time1 = time.time()
    print(f"Processing file {i + 1} of {len(path)}: {file_path}")
    
    # Open the NetCDF file using xarray
    data = xr.open_dataset(file_path)
    
    # Calculate the mean of QGRAUP, QRAIN, and QSNOW along the "bottom_top" dimension
    data['QHYDROMETEORS'] = (data['QGRAUP'] + data['QRAIN'] + data['QSNOW']).mean("bottom_top")
    
    # Drop unnecessary variables
    data1 = data.drop_vars(['QGRAUP', 'QRAIN', 'QSNOW', 'HGT'])
    
    # Convert timestamp to a native Python datetime object
    native_time = time_stamps[i].to_pydatetime()
    
    # Create a DataArray for time with the current timestamp
    time_data_array = xr.DataArray([native_time], dims=['time'], coords={'time': [native_time]})
    
    # Assign the time coordinate to the dataset
    data1 = data1.assign_coords(time=time_data_array)
    
    # Append the processed dataset to the list
    dataset.append(data1)
    
    # Clear memory by deleting variables
    del data, data1, time_data_array
    
    # Collect garbage to release memory
    gc.collect()
    
    # Print a message indicating that memory is cleared for the current file
    print(f"Memory cleared for file {i + 1}: {file_path}")
    
    # Print a message indicating the completion of processing for the current file
    print(f"Processed file {i + 1} of {len(path)}: {file_path}")
    
    # Calculate the time taken to process the current file
    end_time1 = time.time()
    elapsed_time = end_time1 - start_time1
    # Convert seconds to hours, minutes, seconds
    hours, remainder = divmod(elapsed_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    print(f"Processed time for {i+1} of {len(path)}: {file_path} --- {hours} hours, {minutes} minutes, {seconds} seconds")
    print("===================================================================================================================== ")
    
# Concatenate datasets along the time dimension
combine_data = xr.concat(dataset, dim='time')

saved_data_path = "/glade/scratch/noteng/q-hydrometeors.nc"


# Check if "q-hydrometeors.nc" already exists
if os.path.exists(saved_data_path):
    # If it exists, remove it
    print(f"Removing existing NetCDF file: {saved_data_path}")
    os.remove(saved_data_path)

# Save the 'QHYDROMETEORS' variable to a NetCDF file named "q-hydrometeors.nc"
print(f"Saving the combined dataset to {saved_data_path}...")
combine_data['QHYDROMETEORS'].to_netcdf(saved_data_path)

# Print a completion message
print("Processing completed.")

# Calculate the total execution time of the script
end_time = time.time()
elapsed_time = end_time - start_time
# Convert seconds to hours, minutes, seconds
hours, remainder = divmod(elapsed_time, 3600)
minutes, seconds = divmod(remainder, 60)
print("Execution time: ", hours, "hours,", minutes, "minutes,", seconds, "seconds")