#!/bin/bash

# Define the base input directory
base_input_directory="/glade/campaign/ral/wsap/tjuliano/comble/coupled_meso_micro/"

# Define the output directory
output_directory="/glade/scratch/noteng/outfiles1"

# Create the output directory if it doesn't exist
mkdir -p "$output_directory"

# Define a file to keep track of successfully executed files
successful_execution_file="$output_directory/successful_execution.txt"

# Check if the successful execution file exists; if not, create it
if [ ! -f "$successful_execution_file" ]; then
    touch "$successful_execution_file"
fi

# Loop through each RESTART directory
for restart_directory in "$base_input_directory"RESTART_*; do
    # Check if the RESTART directory exists
    if [ -d "$restart_directory" ]; then
        # Loop through each input file matching the pattern
        for input_file in "$restart_directory"/compressed/wrfout_cloud_d02_2020-03-*; do
            # Check if the input file exists
            if [ -f "$input_file" ]; then
                # Extract the file name from the input file path
                input_filename=$(basename "$input_file")
                
                # Create the output file path in the "outfiles" directory with the same name
                output_file="$output_directory/$input_filename"
                
                # Check if the input file has been successfully executed before
                if grep -Fxq "$input_filename" "$successful_execution_file"; then
                    echo "File $input_filename has already been successfully executed. Skipping."
                else
                    # Use ncks to keep only the three variables (QRAIN, QGRAUP, QSNOW)
                    ncks -v "QRAIN,QGRAUP,QSNOW" "$input_file" "$output_file"
                    
                    # Check if ncks was successful
                    if [ $? -eq 0 ]; then
                        echo "Successfully executed $input_filename and saved to $output_file"
                        
                        # Record the successful execution in the tracking file
                        echo "$input_filename" >> "$successful_execution_file"
                    else
                        echo "Error executing $input_filename. Retrying..."
                        
                        # Execute the file again and mark it as successful if successful
                        ncks -v "QRAIN,QGRAUP,QSNOW" "$input_file" "$output_file"
                        if [ $? -eq 0 ]; then
                            echo "Successfully executed $input_filename and saved to $output_file"
                            
                            # Record the successful execution in the tracking file
                            echo "$input_filename" >> "$successful_execution_file"
                        else
                            echo "Error executing $input_filename. Skipping."
                        fi
                    fi
                fi
            else
                echo "Input file $input_file does not exist."
            fi
        done
    else
        echo "RESTART directory $restart_directory does not exist."
    fi
done
