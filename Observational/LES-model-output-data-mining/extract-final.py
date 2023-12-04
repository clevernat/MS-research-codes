import glob
import os
import xarray as xr

# Print import messages
print("Importing Libraries...")

#======================
print("finished importing Libraries")

# Print a message indicating the start of data reading
print("Starting to read data")

# Specify the path to the variables
path_to_variables = sorted(glob.glob("/glade/campaign/ral/wsap/tjuliano/comble/coupled_meso_micro/RESTART_*/compressed/wrfout_cloud_d02_2020-03-*"))

# Specify the path to the latitude and longitude data
path_to_lat_lon = '/glade/campaign/ral/wsap/tjuliano/comble/coupled_meso_micro/RESTART_01/compressed/wrfout_wind_d02_2020-03-13_00_00_00'

# Print a message indicating the end of data reading
print("Done reading data")
#=================

# Print a message indicating the retrieval of height data
print("Getting the height data from the wind path")

# Assuming you've already opened the dataset and have the 'HGT' variable in meters
wind_vars = xr.open_dataset(path_to_lat_lon)
height_meters = wind_vars['HGT']  # 'HGT' is in meters

# Convert height to kilometers
height_kilometers = height_meters / 1000

# Print a message indicating the completion of height conversion
print("Done converting height to kilometers")

# Initialize an empty list to store modified datasets
print("Loading function and performing computations")
datasets_list = []

# Define the output directory
#output_dir = "/glade/scratch/noteng/LES-output/"
output_dir = "/glade/derecho/scratch/noteng/LES-output/"

# Create the output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Loop through the first 10 files in path_to_variables
for i, path in enumerate(path_to_variables[:]):
    # List of variables you want to keep
    variables_to_keep = ['QRAIN', 'QGRAUP', 'QSNOW', 'Times']
    
    # Extract the original file name from the path
    original_filename = os.path.basename(path)
    
    # Define the output file path including the output directory
    output_file = os.path.join(output_dir, original_filename)
    
    # Check if the output file already exists
    if os.path.exists(output_file):
        # Check if the file was successfully executed by checking its size
        if os.path.getsize(output_file) > 0:
            print(f"File '{output_file}' already exists and was successfully executed. Skipping...")
            continue  # Move to the next file
    
    # Open the dataset and perform computations
    try:
        # Open the dataset
        ds = xr.open_dataset(path)
        
        # Remove all variables except the ones in variables_to_keep
        ds = ds.drop_vars([var for var in ds.data_vars if var not in variables_to_keep])
        
        # Add the height_below_1km variable to the dataset
        ds['HGT'] = height_kilometers  # Assuming the dimensions match
        
        # Filter the height variable to include values less than or equal to 1 km 
        ds = ds.where(ds['HGT'] <= 1, drop=True)
        
        # Save the modified dataset to the output file
        ds.to_netcdf(output_file)
        
        # Append the output file name to the list (if needed)
        #datasets_list.append(output_file)
        
        print(f"File '{output_file}' processed and saved.")
    
    except Exception as e:
        print(f"Error processing file '{output_file}': {str(e)}")
        print(f"Retrying file '{output_file}'...")
        continue  # Retry the same file
    
    finally:
        # Close the dataset to release memory
        ds.close()

# Print the list of saved output file names
print(f"Datasets saved to the '{output_dir}' directory.")