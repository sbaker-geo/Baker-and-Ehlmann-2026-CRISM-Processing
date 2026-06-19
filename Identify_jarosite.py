#!/usr/bin/env python
# coding: utf-8

# In[1]:


def detect_jarosite(image_id, ratioed_image, BD2265_image, samples_image, bland_locations):
    
    spectra = Table.read(#path to a .csv file listing the CRISM IF wavelength values, format = 'ascii')
    wl_crism_spectra = spectra['x']
    
    BD2265_list = np.reshape(BD2265_image, BD2265_image.shape[0]*BD2265_image.shape[1])
    
    b3 = ratioed_image[:,:,263]
    b5 = ratioed_image[:,:,266]
    reshaped_b5 = np.reshape(b5, (BD2265_image.shape[0]*BD2265_image.shape[1]))
    reshaped_b3 = np.reshape(b3, (BD2265_image.shape[0]*BD2265_image.shape[1]))
    
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
    

    

    metadata_dict = {}
    metadata_dict['map info'] = '{' + all_map_info[0] + '}'
    metadata_dict['projection info'] = '{' + all_projection_info[0]+ '}'
    metadata_dict['coordinate system string'] = '{' + all_coordinate_system_string[0]+ '}'
    metadata_dict['data ignore value'] = all_data_ignore_value[0]+ '}'
    metadata_dict['default stretch'] = all_default_stretch[0]+ '}'
    
    
    # Gaussian fit
    
    def double_fixed_gauss_func(x, a_1, a_2, c_1, c_2):
        return (-a_1*2.718**(-(((x-4415)**2)/(2*c_1**2)))) + (-a_2*2.718**(-(((x-4535)**2)/(2*c_2**2)))) +1  
    
    double_gauss_data = []
    for i in tqdm.tqdm(range(len(straight_line_cont_rem))):
        if sum(straight_line_cont_rem[i])==33:
            double_gauss_data.append(([0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]))
        elif BD2265_list[i] > 0.003 and reshaped_b5[i] > reshaped_b3[i]:
            try:
                double_gauss_data.append(curve_fit(double_fixed_gauss_func, wn_cropped, straight_line_cont_rem[i], p0=[0.01,0.01,30,30],bounds = ((-np.inf,-np.inf, 0, 0),(np.inf,np.inf, np.inf, np.inf))))
            except:
                double_gauss_data.append(([0,0,0,0],[0,0,0,0],[0,0,0,0], [0,0,0,0]))
        else:
             double_gauss_data.append(([0,0,0,0],[0,0,0,0],[0,0,0,0], [0,0,0,0]))   
    
    popt =[]
    pcov = []
    for i in double_gauss_data:
        popt.append(i[0])
        pcov.append(i[1])              
    a1 = []
    a2 = []
    c1 = []
    c2 = []
    for i in popt:
        a1.append(i[0])
        a2.append(i[1])
        c1.append(i[2])
        c2.append(i[3])
        
    #Save initial depth data as .png and .img
            
    reshaped_a1 = np.reshape(a1, (ratioed_image.shape[0], ratioed_image.shape[1]))
    reshaped_a2 = np.reshape(a2, (ratioed_image.shape[0], ratioed_image.shape[1]))
    reshaped_c1 = np.reshape(c1, (ratioed_image.shape[0], ratioed_image.shape[1]))
    reshaped_c2 = np.reshape(c2, (ratioed_image.shape[0], ratioed_image.shape[1]))
    
    path = #path to save data
    cv2.imwrite(path, reshaped_a1)
    path = #path to save data
    cv2.imwrite(path, reshaped_a2)
    path = #path to save data
    cv2.imwrite(path, reshaped_c1)
    path = #path to save data
    cv2.imwrite(path, reshaped_c2)
    path = #path to save data
    matplotlib.image.imsave(path, reshaped_a1, vmin=-0.03, vmax=0.03, cmap='gray')
    
    output_folder = #path to save data
    base_name = str(image_id) + str('_jar_a1') 
    hyperspectral_data = reshaped_a1
    header_file = os.path.join(output_folder, base_name + '.hdr')
    spectral.envi.save_image(header_file, hyperspectral_data, interleave='bil', dytpe=np.float32, metadata = metadata_dict)
    
    
    #Detect jarosite by comparing to bland standard deviation

    thresholded_a1 = a1.copy()
    
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
            if a1[i] > 3 * mean_ratioed_bland_std:
                thresholded_a1[i] = 1
            else:
                thresholded_a1[i]=0
        return thresholded_a1
           

    def gaussexception(column_number):
        indices = (np.array(np.where(samples_list == column_number)))[0]
        thresholded_gauss_a1 = []
        for i in indices:
            thresholded_a1.append(0)
        return thresholded_a1

    exception_numbers = []
    for number in range(len(cols)):
        try:
            gaussthreshold(number)
        except:
            exception_numbers.append(number)
            
    reshaped_thresholded_a1 = np.reshape(thresholded_a1, (ratioed_image.shape[0], ratioed_image.shape[1]))
    
    column_filtered_jar_classification = []
    for i in range(len(reshaped_thresholded_a1)):
        for j in range(len(reshaped_thresholded_a1[i])):  
            if reshaped_thresholded_a1[i][j]==1 and reshaped_thresholded_a1[i][j+1]==0 and reshaped_thresholded_a1[i][j-1]==0:
                column_filtered_jar_classification.append(0)
            else:
                column_filtered_jar_classification.append(reshaped_thresholded_a1[i][j])
    reshaped_column_filtered_jar_classification = np.reshape(np.array(column_filtered_jar_classification), (ratioed_image.shape[0], ratioed_image.shape[1]))

    path = #path to save data
    cv2.imwrite(path, reshaped_column_filtered_jar_classification)
    output_folder = #path to save data
    base_name = str(image_id) + str('_column_filtered_jar_classification') 
    header_file = os.path.join(output_folder, base_name + '.hdr')
    spectral.envi.save_image(header_file, reshaped_column_filtered_jar_classification, interleave='bil', dytpe=np.float32, metadata = metadata_dict)
    path = #path to save data
    matplotlib.image.imsave(path, reshaped_column_filtered_jar_classification, vmin=0, vmax=1, cmap='gray')
    
    
    column_filtered_thresholded_gauss_a1 = []
    for i in range(len(reshaped_thresholded_a1)):
        for j in range(len(reshaped_thresholded_a1[i])):
            if reshaped_column_filtered_jar_classification[i][j]==1:
                column_filtered_thresholded_gauss_a1.append(reshaped_a1[i][j])
            else:
                column_filtered_thresholded_gauss_a1.append(0)
    reshaped_column_filtered_thresholded_gauss_a1 = np.reshape(np.array(column_filtered_thresholded_gauss_a1), (ratioed_image.shape[0], ratioed_image.shape[1]))
    
    path = #path to save data
    cv2.imwrite(path, reshaped_column_filtered_thresholded_gauss_a1)
    output_folder = #path to save data
    base_name = str(image_id) + str('_column_filtered_thresholded_gauss_a1') 
    header_file = os.path.join(output_folder, base_name + '.hdr')
    spectral.envi.save_image(header_file, reshaped_column_filtered_thresholded_gauss_a1, interleave='bil', dytpe=np.float32, metadata = metadata_dict)
    path = #path to save data
    matplotlib.image.imsave(path, reshaped_column_filtered_thresholded_gauss_a1, vmin=0, vmax=0.03, cmap='gray')
    
    
    path = #path to save data
    np.savetxt(path, a1)
    path = #path to save data
    np.savetxt(path, a2)
    path = #path to save data
    np.savetxt(path, c1)
    path = #path to save data
    np.savetxt(path, c2)
    path = #path to save data
    np.savetxt(path, straight_line_cont_rem)


# In[ ]:




