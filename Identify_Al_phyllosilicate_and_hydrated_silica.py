#!/usr/bin/env python
# coding: utf-8

# In[1]:


def detect_2200_abs(image_id, ratioed_image, LCPINDEX2_image, BD2190_image, BD2250_image, samples_image, bland_locations):
    spectra = Table.read(#path to a .csv file listing the CRISM IF wavelength values, format = 'ascii')
    wl_crism_spectra = spectra['x']
    
    LCPINDEX2_list = np.reshape(LCPINDEX2_image, LCPINDEX2_image.shape[0]*LCPINDEX2_image.shape[1])
    BD2190_list = np.reshape(BD2190_image, LCPINDEX2_image.shape[0]*LCPINDEX2_image.shape[1])
    BD2250_list = np.reshape(BD2250_image, LCPINDEX2_image.shape[0]*LCPINDEX2_image.shape[1])
    
    #Crop to wavelength range
    
    wl_to_use= list(range(240,273))
    cropped_filtered_ratioed_image = ratioed_image[:,:,wl_to_use]
    wn_cropped = (1/wl_crism_spectra[240:273])*10**7
    cropped_filtered_ratioed_image_reshape = np.reshape(np.array(cropped_filtered_ratioed_image), (ratioed_image.shape[0]*ratioed_image.shape[1],33))

    # Straight line continuum remove across this wavelength range
    
    cont_slopes = []
    for i in cropped_filtered_ratioed_image_reshape:
        cont_slopes.append((i[0]-i[-1])/(wn_cropped[0]-wn_cropped[-1]))
    cont_y_int = []
    for i in range(len(cont_slopes)):
        cont_y_int.append(cropped_filtered_ratioed_image_reshape[i][0]-cont_slopes[i]*wn_cropped[0])
    cont_lines = []
    for i in range(len(cont_slopes)):
        cont_lines.append(cont_slopes[i]*np.array(wn_cropped)+cont_y_int[i]) 
    straight_line_cont_rem = np.array(cropped_filtered_ratioed_image_reshape)/np.array(cont_lines)
    
    

    
    # Gaussian fit
        
    
    def single_unfixed_gauss_func(x, a, b, c):
        return (-a*2.718**(-(((x-b)**2)/(2*c**2)))) +1   
    single_gauss_data = []
    for i in tqdm.tqdm(range(len(straight_line_cont_rem))):
        if sum(straight_line_cont_rem[i])==33: #if just straight line don't fit gaussian
            single_gauss_data.append(([0,0,0],[0,0,0],[0,0,0]))
        elif LCPINDEX2_list[i] > 0.027: #if LCP too high don't fit gaussian
            single_gauss_data.append(([0,0,0],[0,0,0],[0,0,0]))
        elif BD2190_list[i] > 0.006 or BD2250_list[i] > 0.006: #if BD high fit gaussian
            try:
                single_gauss_data.append(curve_fit(single_unfixed_gauss_func, wn_cropped, straight_line_cont_rem[i], p0=[0.01, 4500, 30],bounds = ((-np.inf, 4386, 0),(np.inf, 4630, np.inf))))
            except:
                single_gauss_data.append(([0,0,0],[0,0,0],[0,0,0]))
        else:
            single_gauss_data.append(([0,0,0],[0,0,0],[0,0,0]))
            
            
            
    popt =[]
    pcov = []
    for i in single_gauss_data:
        popt.append(i[0])
        pcov.append(i[1])    
    a = []
    b = []
    c = []
    for i in popt:
        a.append(i[0])
        b.append(i[1])
        c.append(i[2])

    #Save initial depth data as .png and .img
    
    reshaped_a = np.reshape(a, (ratioed_image.shape[0],ratioed_image.shape[1]))
    path = #path to save data
    matplotlib.image.imsave(path, reshaped_a, vmin=-0.03, vmax=0.03, cmap='gray')
    
    metadata_dict = {}
    metadata_dict['map info'] = '{' + all_map_info[0] + '}'
    metadata_dict['projection info'] = '{' + all_projection_info[0]+ '}'
    metadata_dict['coordinate system string'] = '{' + all_coordinate_system_string[0]+ '}'
    metadata_dict['data ignore value'] = all_data_ignore_value[0]+ '}'
    metadata_dict['default stretch'] = all_default_stretch[0]+ '}'
    
    output_folder = #path to save data'
    base_name = str(image_id) + str('_2200_abs_a') 
    hyperspectral_data = reshaped_a
    header_file = os.path.join(output_folder, base_name + '.hdr')
    spectral.envi.save_image(header_file, hyperspectral_data, interleave='bil', dytpe=np.float32, metadata = metadata_dict)
    
    reshaped_b = np.reshape(b, (ratioed_image.shape[0],ratioed_image.shape[1],1))
    
    #Detect 2200 abs by comparing to bland standard deviation
    
    thresholded_a = a.copy()
    
    samples_list = np.reshape(samples_image, samples_image.shape[0]*samples_image.shape[1])
    max_samples = int(max(samples_list)[0])
    cols = np.linspace(0,max_samples,max_samples - 1) 
    location_of_blandest_reshaped = np.reshape(bland_locations, (ratioed_image.shape[0]*ratioed_image.shape[1]))
    
    def gaussthreshold(column_number):
        indices = (np.array(np.where(samples_list == column_number)))[0]
        bland_std_in_col = []
        for i in indices:
            if location_of_blandest_reshaped[i] == 1:
                bland_std_in_col.append(np.std(np.array(straight_line_cont_rem[i])))
        mean_ratioed_bland_std = np.mean(bland_std_in_col)
        for i in indices:
            if a[i] > 3 * mean_ratioed_bland_std:
                thresholded_a[i] = 1
            else:
                thresholded_a[i]=0
        return thresholded_a
   

    def gaussexception(column_number):
        indices = (np.array(np.where(samples_list == column_number)))[0]
        thresholded_gauss_a = []
        for i in indices:
            thresholded_a.append(0)
        return thresholded_a

    exception_numbers = []
    for number in range(len(cols)):
        try:
            gaussthreshold(number)
        except:
            exception_numbers.append(number)
    
    reshaped_thresholded_a = np.reshape(thresholded_a, (ratioed_image.shape[0],ratioed_image.shape[1]))
    
    column_filtered_2200_classification = []
    for i in range(len(reshaped_thresholded_a)):
        for j in range(len(reshaped_thresholded_a[i])):
            if reshaped_thresholded_a[i][j]==1 and reshaped_thresholded_a[i][j+1]==0 and reshaped_thresholded_a[i][j-1]==0:
                column_filtered_2200_classification.append(0)
            elif reshaped_b[i][j] < 4444.4:
                column_filtered_2200_classification.append(0)
            else:
                column_filtered_2200_classification.append(reshaped_thresholded_a[i][j])
    reshaped_column_filtered_2200_classification = np.reshape(np.array(column_filtered_2200_classification), (ratioed_image.shape[0],ratioed_image.shape[1]))

    
    path = #path to save data
    cv2.imwrite(path, reshaped_column_filtered_2200_classification)
    output_folder = #path to save data
    base_name = str(image_id) + str('_column_filtered_2200_classification') 
    header_file = os.path.join(output_folder, base_name + '.hdr')
    spectral.envi.save_image(header_file, reshaped_column_filtered_2200_classification, interleave='bil', dytpe=np.float32, metadata = metadata_dict)
    path = #path to save data
    matplotlib.image.imsave(path, reshaped_column_filtered_2200_classification, vmin=0, vmax=1, cmap='gray')
    
    
    column_filtered_thresholded_gauss_a = []
    for i in range(len(reshaped_thresholded_a)):
        for j in range(len(reshaped_thresholded_a[i])):
            if reshaped_column_filtered_2200_classification[i][j]==1:
                column_filtered_thresholded_gauss_a.append(reshaped_a[i][j])
            else:
                column_filtered_thresholded_gauss_a.append(0)
    reshaped_column_filtered_thresholded_gauss_a = np.reshape(np.array(column_filtered_thresholded_gauss_a), (ratioed_image.shape[0],ratioed_image.shape[1]))
    path = #path to save data
    cv2.imwrite(path, reshaped_column_filtered_thresholded_gauss_a)
    output_folder = #path to save data
    base_name = str(image_id) + str('_column_filtered_thresholded_gauss_a') 
    header_file = os.path.join(output_folder, base_name + '.hdr')
    spectral.envi.save_image(header_file, reshaped_column_filtered_thresholded_gauss_a, interleave='bil', dytpe=np.float32, metadata = metadata_dict)
    path = #path to save data
    matplotlib.image.imsave(path, reshaped_column_filtered_thresholded_gauss_a, vmin=0, vmax=0.03, cmap='gray')
    
    
    
    path = #path to save data
    np.savetxt(path, a)
    path = #path to save data
    np.savetxt(path, b)
    path = #path to save data
    np.savetxt(path, c)
    path = #path to save data
    np.savetxt(path, straight_line_cont_rem)
    


# In[ ]:




