#!/usr/bin/env python
# coding: utf-8

# In[1]:


def ratio_image(image_id, masked_image, samples, BD530_2, BD860_2, BD1400, OLINDEX3, LCPINDEX2, HCPINDEX2,
                BD1750_2, BD1900_2, BD2190, BD2250, BD2265, BD2290, BD2355):
    spectra = Table.read(#path to a .csv file listing the CRISM IF wavelength values, format = 'ascii')
    wl_crism_spectra = spectra['wavelength']
    
    # Continuum remove and calculating the area under the curve
    
    bands_to_use= list(range(203,286))
    whole_image = masked_image[:,:,bands_to_use]
    wn_cropped = (1/wl_crism_spectra[203:286])*10**7
    whole_image_reshape = np.reshape(np.array(whole_image), (whole_image.shape[0]*whole_image.shape[1],83))
    

    def continuum_removal(points, show=False):
        x, y = points.T
        augmented = np.concatenate([points, [(x[0], np.min(y)-1), (x[-1], np.min(y)-1)]], axis=0)
        hull = ConvexHull(augmented)
        continuum_points = points[np.sort([v for v in hull.vertices if v < len(points)])]
        continuum_function = interp1d(*continuum_points.T)
        yprime = y / continuum_function(x)
        if show:
            fig, axes = plt.subplots(2, 1, sharex=True)
            axes[0].plot(x, y, label='Data')
            axes[0].plot(*continuum_points.T, label='Continuum')
            axes[0].legend()
            axes[1].plot(x, yprime, label='Data / Continuum')
            axes[1].legend()
        return np.c_[x, yprime],continuum_points, yprime


    points = []
    for i in range(len(whole_image_reshape)):
        points.append(np.c_[wn_cropped, whole_image_reshape[i]])  
    points=np.array(points)               
    contrem_data = []
    for i in range(len(points)):
        contrem_data.append((continuum_removal(points[i], show = False)))
                
    new_points = [] 
    continuum_points = []
    yprime = [] 

    for i in contrem_data:
        new_points.append(i[0])
        continuum_points.append(i[1])
        yprime.append(i[2])

        
    area_cont = (wn_cropped[0]-wn_cropped[-1])*1
    area_spec = []
    for i in range(len(yprime)):
        area_spec.append(abs(np.trapz(yprime[i], x=wn_cropped)))
    
    area_abs = []
    for i in range(len(area_spec)):
        area_abs.append(area_cont - area_spec[i])
    
    #Make each parameter a 1D list
        
    BD530_2_list = BD530_2.reshape(BD530_2.shape[0]*BD530_2.shape[1])
    BD860_2_list = BD860_2.reshape(BD860_2.shape[0]*BD530_2.shape[1])
    BD1400_list = BD1400.reshape(BD1400.shape[0]*BD530_2.shape[1])
    OLINDEX3_list = OLINDEX3.reshape(OLINDEX3.shape[0]*BD530_2.shape[1])
    LCPINDEX2_list = LCPINDEX2.reshape(LCPINDEX2.shape[0]*BD530_2.shape[1])
    HCPINDEX2_list = HCPINDEX2.reshape(HCPINDEX2.shape[0]*BD530_2.shape[1])
    BD1750_2_list = BD1750_2.reshape(BD1750_2.shape[0]*BD530_2.shape[1])
    BD1900_2_list = BD1900_2.reshape(BD1900_2.shape[0]*BD530_2.shape[1])
    BD2190_list = BD2190.reshape(BD2190.shape[0]*BD530_2.shape[1])
    BD2250_list = BD2250.reshape(BD2250.shape[0]*BD530_2.shape[1])
    BD2265_list = BD2265.reshape(BD2265.shape[0]*BD530_2.shape[1])
    BD2290_list = BD2290.reshape(BD2290.shape[0]*BD530_2.shape[1])
    BD2355_list = BD2355.reshape(BD2355.shape[0]*BD2355.shape[1])
    
    
    
    samples_list = samples.reshape(samples.shape[0]*samples.shape[1])
    max_samples = int(max(samples_list)[0])
    reshaped_masked_image = np.reshape(masked_image, (masked_image.shape[0]*masked_image.shape[1], 489))
    ratioed_image = reshaped_masked_image.copy()
    cols = np.linspace(0,max_samples,max_samples - 1) 
    

    def findBlandest(column_number):
        indices = (np.array(np.where(samples_list == column_number)))[0]
        BD530_2_in_column = []
        BD860_2_in_column = []
        BD1400_in_column = []
        OLINDEX3_in_column = []
        LCPINDEX2_in_column = []
        HCPINDEX2_in_column = []
        BD1750_2_in_column = []
        BD1900_2_in_column = []
        BD2190_in_column = []
        BD2250_in_column = []
        BD2265_in_column = []
        BD2290_in_column = []
        BD2355_in_column = []
        for i in indices:
            BD530_2_in_column.append((BD530_2_list[i]))
            BD860_2_in_column.append(BD860_2_list[i])
            BD1400_in_column.append(BD1400_list[i])
            OLINDEX3_in_column.append(OLINDEX3_list[i])
            LCPINDEX2_in_column.append(LCPINDEX2_list[i])
            HCPINDEX2_in_column.append(HCPINDEX2_list[i])
            BD1750_2_in_column.append(BD1750_2_list[i])
            BD1900_2_in_column.append(BD1900_2_list[i])
            BD2190_in_column.append(BD2190_list[i])
            BD2250_in_column.append(BD2250_list[i])
            BD2265_in_column.append(BD2265_list[i])
            BD2290_in_column.append(BD2290_list[i])
            BD2355_in_column.append(BD2355_list[i])
        
        BD530_2_in_column_real = []
        BD860_2_in_column_real = []
        BD1400_in_column_real = []
        OLINDEX3_in_column_real = []
        LCPINDEX2_in_column_real = []
        HCPINDEX2_in_column_real = []
        BD1750_2_in_column_real = []
        BD1900_2_in_column_real = []
        BD2190_in_column_real = []
        BD2250_in_column_real = []
        BD2265_in_column_real = []
        BD2290_in_column_real = []
        BD2355_in_column_real = []
        for i in range(len(BD530_2_in_column)):
            if BD530_2_in_column[i] != '--':
                BD530_2_in_column_real.append(float(BD530_2_in_column[i]))
                BD860_2_in_column_real.append(float(BD860_2_in_column[i]))
                BD1400_in_column_real.append(float(BD1400_in_column[i]))
                OLINDEX3_in_column_real.append(float(OLINDEX3_in_column[i]))
                LCPINDEX2_in_column_real.append(float(LCPINDEX2_in_column[i]))
                HCPINDEX2_in_column_real.append(float(HCPINDEX2_in_column[i]))
                BD1750_2_in_column_real.append(float(BD1750_2_in_column[i]))
                BD1900_2_in_column_real.append(float(BD1900_2_in_column[i]))
                BD2190_in_column_real.append(float(BD2190_in_column[i]))
                BD2250_in_column_real.append(float(BD2250_in_column[i]))
                BD2265_in_column_real.append(float(BD2265_in_column[i]))
                BD2290_in_column_real.append(float(BD2290_in_column[i]))
                BD2355_in_column_real.append(float(BD2355_in_column[i]))

        index_of_blandest = []
        for i in range(len(BD530_2_in_column)):
            if (BD530_2_in_column[i] < (statistics.median(BD530_2_in_column_real) + 1.5*(np.std(BD530_2_in_column_real))))and ((BD860_2_in_column[i] < statistics.median(BD860_2_in_column_real) + 1.5*(np.std(BD860_2_in_column_real))))and ((BD860_2_in_column[i] < statistics.median(BD860_2_in_column_real) + 1.5*(np.std(BD860_2_in_column_real))))and ((BD1400_in_column[i] < statistics.median(BD1400_in_column_real) + 1.5*(np.std(BD1400_in_column_real))))and ((OLINDEX3_in_column[i] < statistics.median(OLINDEX3_in_column_real) + 1.5*(np.std(OLINDEX3_in_column_real))))and ((LCPINDEX2_in_column[i] < statistics.median(LCPINDEX2_in_column_real) + 1.5*(np.std(LCPINDEX2_in_column_real))))and ((HCPINDEX2_in_column[i] < statistics.median(HCPINDEX2_in_column_real) + 1.5*(np.std(HCPINDEX2_in_column_real))))and ((BD1750_2_in_column[i] < statistics.median(BD1750_2_in_column_real) + 1.5*(np.std(BD1750_2_in_column_real))))and ((BD1900_2_in_column[i] < statistics.median(BD1900_2_in_column_real) + 1.5*(np.std(BD1900_2_in_column_real))))and ((BD2190_in_column[i] < statistics.median(BD2190_in_column_real) + 1.5*(np.std(BD2190_in_column_real))))and ((BD2250_in_column[i] < statistics.median(BD2250_in_column_real) + 1.5*(np.std(BD2250_in_column_real))))and ((BD2265_in_column[i] < statistics.median(BD2265_in_column_real) + 1.5*(np.std(BD2265_in_column_real))))and ((BD2290_in_column[i] < statistics.median(BD2290_in_column_real) + 1.5*(np.std(BD2290_in_column_real))))and ((BD2355_in_column[i] < statistics.median(BD2355_in_column_real) + 1.5*(np.std(BD2355_in_column_real)))) and (BD2290_in_column[i] < 0.003):
                index_of_blandest.append(indices[i])
        specs_in_column = []
        for i in index_of_blandest:
            if area_abs[i] !=0:
                specs_in_column.append(area_abs[i])
        sorted_specs_in_column = np.sort(specs_in_column)
        blandest_areas = sorted_specs_in_column[0:9]
        index_of_blandest_unmerged = []
        for i in blandest_areas:
            index_of_blandest_unmerged.append(list(np.where(area_abs == i)[0]))
        new_index_of_blandest = list(itertools.chain.from_iterable(index_of_blandest_unmerged))
        blandest_specs = []
        for i in new_index_of_blandest:
            blandest_specs.append(np.array(reshaped_masked_image[i]))
        mean_blandest = sum(blandest_specs)/len(blandest_specs)
        for i in indices:
            ratioed_image[i]=(reshaped_masked_image[i]/mean_blandest)
        return ratioed_image

    
    def exception1(column_number):
        indices = (np.array(np.where(samples_list == column_number)))[0]           
        blandest_specs = []
        for i in indices:
            blandest_specs.append(np.array(reshaped_masked_image[i]))
        filtered_blandest_specs = []
        for i in blandest_specs:
            if 65535 not in i:
                filtered_blandest_specs.append(i)
        mean_blandest = sum(filtered_blandest_specs)/len(filtered_blandest_specs)
        for i in indices:
                ratioed_image[i]=(reshaped_masked_image[i]/mean_blandest)
        return ratioed_image
    
    
    def exception2(column_number):
        indices = (np.array(np.where(samples_list == column_number)))[0]
        for i in indices:
            ratioed_image[i]=(np.linspace(0,0,489))
        return ratioed_image

        
     
    for number in range(len(cols)):
        try:
            findBlandest(number)
        except:
            try:
                exception1(number)
            except:
                exception2(number)
     

        
    reshaped_ratioed_image = []
    for i in ratioed_image:
        reshaped_ratioed_image.append(np.reshape(i, 489))
        
    re_reshaped_ratioed_image = np.reshape(reshaped_ratioed_image, (whole_image.shape[0],whole_image.shape[1],489))
    
    path = #path to save data
    np.savetxt(path, reshaped_ratioed_image)
        
    
    SG_filtered_ratioed = []
    for i in range(len(reshaped_ratioed_image)):
        try:
            SG_filtered_ratioed.append(scipy.signal.savgol_filter(reshaped_ratioed_image[i], 5,2))
        except:
            SG_filtered_ratioed.append(list([0]*489))    
    reshaped_SG_filtered_ratioed = np.reshape(SG_filtered_ratioed, (whole_image.shape[0],whole_image.shape[1],489))
    
    
    metadata_dict_wavelength = {}
    metadata_dict_wavelength['map info'] = '{' + all_map_info[0] + '}'
    metadata_dict_wavelength['projection info'] = '{' + all_projection_info[0]+ '}'
    metadata_dict_wavelength['coordinate system string'] = '{' + all_coordinate_system_string[0]+ '}'
    metadata_dict_wavelength['default bands'] = all_default_bands[0]+ '}'
    metadata_dict_wavelength['wavelength units'] = all_wavelength_units[0]+ '}'
    metadata_dict_wavelength['data ignore value'] = all_data_ignore_value[0]+ '}'
    metadata_dict_wavelength['default stretch'] = all_default_stretch[0]+ '}'
    metadata_dict_wavelength['wavelength'] = '{' + all_wavelength[0]+ '}'
    metadata_dict_wavelength['fwhm'] = '{' + all_fwhm[0]+ '}'
    metadata_dict_wavelength['bbl'] = '{' + all_bbl[0] + '}'
    
    
    output_folder = #path to save data
    base_name = str(image_id) + str('_ratioed_SG_image') 
    hyperspectral_data = reshaped_SG_filtered_ratioed
    header_file = os.path.join(output_folder, base_name + '.hdr')
    spectral.envi.save_image(header_file, hyperspectral_data, interleave='bil', dytpe=np.float32, metadata = metadata_dict_wavelength)
    
    return hyperspectral_data
    
    


# In[ ]:




