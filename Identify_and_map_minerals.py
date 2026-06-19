#!/usr/bin/env python
# coding: utf-8

# In[2]:


#Importing functions

get_ipython().run_line_magic('run', 'Identify_FeMg_phyllosilicate.ipynb')
get_ipython().run_line_magic('run', 'Identify_Fe_oxide.ipynb')
get_ipython().run_line_magic('run', 'Identify_jarosite.ipynb')
get_ipython().run_line_magic('run', 'Identify_Al_phyllosilicate_and_hydrated_silica.ipynb')
get_ipython().run_line_magic('run', 'Create_ratioed_images.ipynb')
get_ipython().run_line_magic('run', 'Identify_location_of_blandest_pixels.ipynb')


# In[3]:


# Importing usual functions

import numpy as np
import matplotlib.pyplot as plt
import pathlib as pathlib
from pathlib import Path
import os
from os import getcwd, path, chdir, makedirs, listdir
from time import strftime
from csv import writer
import spectral
from spectral import open_image, calc_stats, noise_from_diffs, mnf, envi
from arsf_envi_reader import envi_header as eh
import statistics
import PIL
from PIL import Image
import cv2
from astropy.table import Table
import scipy
from scipy.optimize import curve_fit
import itertools
from scipy.spatial import ConvexHull
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import tqdm
import statistics
import matplotlib.image
import rasterio
from skimage import morphology, measure
from scipy import stats
import glob
import gc
import math


# In[4]:


# Read in all files names in a folder of CRISM images

folder_path = #path to your folder

file_names = []
for item in os.listdir(folder_path):
    item_path = os.path.join(folder_path, item)
    if os.path.isfile(item_path) and 'if' in item and 'hdr' in item:
        file_names.append(item)

print(file_names)


# In[5]:


# May need to change the indices to crop the image ID

image_names = []

for i in file_names:
    if i[6] == 0:
        image_names.append(i[7:11])
    else:
        image_names.append(i[7:11])

print(image_names)   


# In[6]:


for image in image_names:
    
    image_to_use= str(image)
    
    # %% SET VARIABLES
    dataDir = #Path to your folder of CRISM images
    subset = ['waves', 1000, 2600]

    # %% DEFINE FUNCTIONS AND CLASSES
    class Directory:
        def __init__(self, dataDir):
            self.data = Path(dataDir)
            self.code = getcwd()
            
        
    def get_header_names(directory):
        return([file for file in listdir(directory)                 if 'if' in file and image in file and file.endswith('.hdr') or file.endswith('.HDR')])      
                

    def get_header(directory, headerName, subset):
        chdir(directory)
        header = eh.read_hdr_file(headerName)
        wavelengths = [float(w) for w in header['wavelength'].split(',')]
        header.wavelength = np.array(wavelengths)
        if subset[0].lower() == 'bands':
            header['bandList'] = list(range(subset[1]-1, subset[2]))
        elif subset[0].lower() == 'waves':
            min = np.argmin(np.abs(np.array(header.wavelength-subset[1])))
            max = np.argmin(np.abs(np.array(header.wavelength-subset[2])))
            header['bandList'] = list(range(min, max))
        elif subset[0].lower() == 'all':
            header['bandList'] = list(range(0, len(wavelengths)))
        else:
            raise Exception("Subset must be 'bands', 'waves', or 'all'")
        return(header)


    headerNames = get_header_names(dataDir)

    with open('nf_images.csv_2','w') as f:
        for row in headerNames:
            for x in row:
                f.write(str(x))
            f.write('\n')
        

    # %% MAIN SCRIPT
    dir = Directory(dataDir)
    headerNames = get_header_names(dir.data)

    headerNames.sort()

    all_masked_image = []

    for headerName in headerNames:
        header = get_header(dir.data, headerName, subset)
        image = open_image(headerName).load()
        all_masked_image.append(np.ma.MaskedArray(data=image, mask=(image == 65535.0)))
        
  
    #**************************************************************************************************
    # Get metadata

    def get_metadata(directory, headerName, subset):
        chdir(directory)
        header = eh.read_hdr_file(headerName)
        map_info = header['map info']
        projection_info = header['projection info']
        coordinate_system_string = header['coordinate system string']
        default_bands = header['default bands']
        wavelength_units = header['wavelength units']
        data_ignore_value = header['data ignore value']
        default_stretch = header['default stretch']
        wavelength = header['wavelength']
        fwhm = header['fwhm']
        bbl = header['bbl']
        return map_info, projection_info, coordinate_system_string, default_bands, wavelength_units, data_ignore_value, default_stretch, wavelength, fwhm, bbl

    all_map_info = []
    all_projection_info = []
    all_coordinate_system_string = []
    all_default_bands = []
    all_wavelength_units = []
    all_data_ignore_value = []
    all_default_stretch = []
    all_wavelength = []
    all_fwhm = []
    all_bbl = []

    for headerName in headerNames:
        header = get_header(dir.data, headerName, subset)
        map_info, projection_info, coordinate_system_string, default_bands, wavelength_units, data_ignore_value, default_stretch, wavelength, fwhm, bbl = get_metadata(dir.data, headerName, subset)
        all_map_info.append(map_info)
        all_projection_info.append(projection_info)
        all_coordinate_system_string.append(coordinate_system_string)
        all_default_bands.append(default_bands)
        all_wavelength_units.append(wavelength_units)
        all_data_ignore_value.append(data_ignore_value)
        all_default_stretch.append(default_stretch)
        all_wavelength.append(wavelength)
        all_fwhm.append(fwhm)
        all_bbl.append(bbl)
    
  
    #**************************************************************************************************
    # Get su data
    
    
    def get_header_names(directory):
        return([file for file in listdir(directory)                 if 'su' in file and image_to_use in file and file.endswith('.hdr') or file.endswith('.HDR')])


    def get_header(directory, headerName, subset):
        chdir(directory)
        header = eh.read_hdr_file(headerName)
        wavelengths = [float(w) for w in header['band names'].split(',')]
        header.wavelength = np.array(wavelengths)
        if subset[0].lower() == 'bands':
            header['bandList'] = list(range(subset[1]-1, subset[2]))
        elif subset[0].lower() == 'waves':
            min = np.argmin(np.abs(np.array(header.wavelength-subset[1])))
            max = np.argmin(np.abs(np.array(header.wavelength-subset[2])))
            header['bandList'] = list(range(min, max))
        elif subset[0].lower() == 'all':
            header['bandList'] = list(range(0, len(wavelengths)))
        else:
            raise Exception("Subset must be 'bands', 'waves', or 'all'")
        return(header)
    headerNames = get_header_names(dataDir)
    dir = Directory(dataDir)
    headerNames = get_header_names(dir.data)
    headerNames.sort()
    all_masked_image_2 = []
    for headerName in headerNames: 
        image = open_image(headerName).load()
        masked_image_2 = np.ma.MaskedArray(data=image, mask=(image == 65535.0))
        all_masked_image_2.append(masked_image_2)
    
    
    all_BD530_2 = []
    all_BD860_2 = []
    all_BD1400 = []
    all_OLINDEX3 = []
    all_LCPINDEX2 = []
    all_HCPINDEX2 = []
    all_BD1750_2 = []
    all_BD1900_2 = []
    all_BD2190 = []
    all_BD2250 = []
    all_BD2265 = []
    all_BD2290 = []
    all_BD2355 = []

    for masked_image_2 in all_masked_image_2:
        all_BD530_2.append(masked_image_2[:,:,2])
        all_BD860_2.append(masked_image_2[:,:,6])
        all_BD1400.append(masked_image_2[:,:,20])
        all_OLINDEX3.append(masked_image_2[:,:,13])
        all_LCPINDEX2.append(masked_image_2[:,:,16])
        all_HCPINDEX2.append(masked_image_2[:,:,17])
        all_BD1750_2.append(masked_image_2[:,:,24])
        all_BD1900_2.append(masked_image_2[:,:,25])
        all_BD2190.append(masked_image_2[:,:,30])
        all_BD2250.append(masked_image_2[:,:,35])
        all_BD2265.append(masked_image_2[:,:,37])
        all_BD2290.append(masked_image_2[:,:,38])
        all_BD2355.append(masked_image_2[:,:,40])
    
    
    #**************************************************************************************************
    # Get IN data 
    
    
    def get_header_names(directory):
        return([file for file in listdir(directory)             if 'in' in file and image_to_use in file and file.endswith('.hdr') or file.endswith('.HDR')])
    def get_header(directory, headerName, subset):
        chdir(directory)
        header = eh.read_hdr_file(headerName)
        wavelengths = [float(w) for w in header['band names'].split(',')]
        header.wavelength = np.array(wavelengths)
        if subset[0].lower() == 'bands':
            header['bandList'] = list(range(subset[1]-1, subset[2]))
        elif subset[0].lower() == 'waves':
            min = np.argmin(np.abs(np.array(header.wavelength-subset[1])))
            max = np.argmin(np.abs(np.array(header.wavelength-subset[2])))
            header['bandList'] = list(range(min, max))
        elif subset[0].lower() == 'all':
            header['bandList'] = list(range(0, len(wavelengths)))
        else:
            raise Exception("Subset must be 'bands', 'waves', or 'all'")
        return(header)
    headerNames = get_header_names(dataDir)
    dir = Directory(dataDir)
    headerNames = get_header_names(dir.data)
    headerNames.sort()
    all_samples = []
    for headerName in headerNames:
        image = open_image(headerName).load()
        masked_image_2 = np.ma.MaskedArray(data=image, mask=(image == 65535.0))
        all_samples.append(masked_image_2[:,:,3])
    
    
    image = image_to_use
    
    #**************************************************************************************************
    # Run functions
    
    ratioed_spectra_image = ratio_image(image_to_use, all_masked_image[0], all_samples[0], all_BD530_2[0],
                                    all_BD860_2[0],all_BD1400[0], all_OLINDEX3[0], all_LCPINDEX2[0],
                                    all_HCPINDEX2[0], all_BD1750_2[0], all_BD1900_2[0], all_BD2190[0],
                                    all_BD2250[0], all_BD2265[0], all_BD2290[0], all_BD2355[0])
    
    location_of_bland_image = location_of_bland(image_to_use, all_masked_image[0], all_samples[0], all_BD530_2[0],
                                    all_BD860_2[0],all_BD1400[0], all_OLINDEX3[0], all_LCPINDEX2[0],
                                    all_HCPINDEX2[0], all_BD1750_2[0], all_BD1900_2[0], all_BD2190[0],
                                    all_BD2250[0], all_BD2265[0], all_BD2290[0], all_BD2355[0])
    
    detect_2200_abs(image_to_use, ratioed_spectra_image, all_LCPINDEX2[0], all_BD2190[0], all_BD2250[0],
                    all_samples[0],location_of_bland_image)
    
    detect_FeMg_phyllo(image_to_use, ratioed_spectra_image, all_BD2290[0], all_samples[0], location_of_bland_image)
    
    detect_jarosite(image_to_use, ratioed_spectra_image, all_BD2265[0], all_samples[0], location_of_bland_image)
    
    detect_Fe_oxide(image_to_use, ratioed_spectra_image, all_BD530_2[0], all_samples[0], location_of_bland_image)

    #Clear variables before running the loop again
    
    del all_masked_image
    del all_masked_image_2
    del ratioed_spectra_image
   


# In[ ]:




