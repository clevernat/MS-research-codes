import numpy as np

### Function to find position in 1d array closest to a value ###
def geo_idx(dd, dd_array):
    geo_idx = (np.abs(dd_array - dd)).argmin()
    return geo_idx
