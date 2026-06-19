#!/usr/bin/env python
# coding: utf-8

# In[1]:


def detect_Fe_oxide(image_id, ratioed_image, BD530_image, samples_image, bland_locations):
    
    spectra = Table.read(#path to a .csv file listing the CRISM IF wavelength values, format = 'ascii')
    wl_crism_spectra = spectra['x']
    
    BD530_2_list = np.reshape(BD530_image, BD530_image.shape[0]*BD530_image.shape[1])
    

    #Crop to wavelength range
       
    wl_to_use= list(range(1,30))
    cropped_filtered_ratioed_image = ratioed_image[:,:,wl_to_use]
    wn_cropped_530 = (1/wl_crism_spectra[1:30])*10**7
    cropped_filtered_ratioed_image_reshape_530 = np.reshape(np.array(cropped_filtered_ratioed_image), (ratioed_image.shape[0]*ratioed_image.shape[1],29))
    
    # Straight line continuum remove across this wavelength range
    
    cont_slopes_530 = []
    for i in cropped_filtered_ratioed_image_reshape_530:
        cont_slopes_530.append((i[0]-i[-1])/(wn_cropped_530[0]-wn_cropped_530[-1]))
    cont_y_int_530 = []
    for i in range(len(cont_slopes_530)):
        cont_y_int_530.append(cropped_filtered_ratioed_image_reshape_530[i][0]-cont_slopes_530[i]*wn_cropped_530[0])
    cont_lines_530 = []
    for i in range(len(cont_slopes_530)):
        cont_lines_530.append(cont_slopes_530[i]*np.array(wn_cropped_530)+cont_y_int_530[i]) 
    straight_line_cont_rem_530 = np.array(cropped_filtered_ratioed_image_reshape_530)/np.array(cont_lines_530)  
    
    
    # Gaussian fit
    
    def single_unfixed_gauss_func(x, a, b, c):
        return (-a*2.718**(-(((x-b)**2)/(2*c**2)))) +1
    
    single_gauss_data_530 = []

    for i in tqdm.tqdm(range(len(straight_line_cont_rem_530))):
        if sum(straight_line_cont_rem_530[i])==29:
            single_gauss_data_530.append(([0,0,0],[0,0,0],[0,0,0]))
        elif BD530_2_list[i] > 0.16:
            try:
                single_gauss_data_530.append(curve_fit(single_unfixed_gauss_func, wn_cropped_530, straight_line_cont_rem_530[i], p0=[0.03, 20000, 1000],bounds = ((-np.inf, 16667, 0),(np.inf, 25000, np.inf))))
            except:
                single_gauss_data_530.append(([0,0,0],[0,0,0],[0,0,0]))
        else:
            single_gauss_data_530.append(([0,0,0],[0,0,0],[0,0,0]))
            
            
            
    popt_530 =[]
    pcov_530 = []
    for i in single_gauss_data_530:
        popt_530.append(i[0])
        pcov_530.append(i[1])    
    a_530 = []
    b_530 = []
    c_530 = []
    for i in popt_530:
        a_530.append(i[0])
        b_530.append(i[1])
        c_530.append(i[2])
    
    #Save initial depth data as .png and .img
    
    metadata_dict = {}
    metadata_dict['map info'] = '{' + all_map_info[0] + '}'
    metadata_dict['projection info'] = '{' + all_projection_info[0]+ '}'
    metadata_dict['coordinate system string'] = '{' + all_coordinate_system_string[0]+ '}'
    metadata_dict['data ignore value'] = all_data_ignore_value[0]+ '}'
    metadata_dict['default stretch'] = all_default_stretch[0]+ '}'
    
    reshaped_a_530 = np.reshape(a_530, (ratioed_image.shape[0],ratioed_image.shape[1]))
    path = #path to save data
    cv2.imwrite(path, reshaped_a_530)
    path = #path to save data
    matplotlib.image.imsave(path, reshaped_a_530, vmin=-0.03, vmax=0.03, cmap='gray')
    output_folder = '/export/data1/sbaker/LM_data/LM_Fe_oxide/'
    base_name = str(image_id) + str('_Fe_oxide_a_530') 
    hyperspectral_data = reshaped_a_530
    header_file = os.path.join(output_folder, base_name + '.hdr')
    spectral.envi.save_image(header_file, hyperspectral_data, interleave='bil', dytpe=np.float32, metadata = metadata_dict)
    
    
    reshaped_b_530 = np.reshape(b_530, (ratioed_image.shape[0],ratioed_image.shape[1],1))
    
    #Detect 2200 abs by comparing to bland standard deviation

    thresholded_a_530 = a_530.copy()
    
    samples_list = np.reshape(samples_image, samples_image.shape[0]*samples_image.shape[1])
    max_samples = int(max(samples_list)[0])
    cols = np.linspace(0,max_samples,max_samples - 1) 
    location_of_blandest_reshaped = np.reshape(bland_locations, (ratioed_image.shape[0]*ratioed_image.shape[1]))
    
    def gaussthreshold_530(column_number):
        indices = (np.array(np.where(samples_list == column_number)))[0]
        bland_std_in_col_530 = []
        for i in indices:
            if location_of_blandest_reshaped[i] == 1:
                bland_std_in_col_530.append(np.std(np.array(straight_line_cont_rem_530[i])))
        mean_ratioed_bland_std_530 = np.mean(bland_std_in_col_530)
        for i in indices:
            if a_530[i] > 3 * mean_ratioed_bland_std_530:
                thresholded_a_530[i] = 1
            else:
                thresholded_a_530[i]=0
        return thresholded_a_530 

       

    def gaussexception(column_number):
        indices = (np.array(np.where(samples_list == column_number)))[0]
        thresholded_gauss_a = []
        for i in indices:
            thresholded_a.append(0)
        return thresholded_a

    exception_numbers = []
    for number in range(len(cols)):
        try:
            gaussthreshold_530(number)
        except:
            exception_numbers.append(number)

    reshaped_thresholded_a_530 = np.reshape(thresholded_a_530, (ratioed_image.shape[0],ratioed_image.shape[1]))
    
    
    wl_to_use= list(range(0,75))
    cropped_filtered_ratioed_image = ratioed_image[:,:,wl_to_use]
    wn_cropped = (1/wl_crism_spectra[0:75])*10**7
    cropped_filtered_ratioed_image_reshape = np.reshape(np.array(cropped_filtered_ratioed_image), (ratioed_image.shape[0]*ratioed_image.shape[1],75))


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
    
    BD625 = []
    for i in range(len(straight_line_cont_rem)):
        BD625.append(1 - straight_line_cont_rem[i][29]/(0.667*straight_line_cont_rem[i][0] + 0.333*straight_line_cont_rem[i][74]))
    
    BD625_reshaped = np.reshape(BD625, ((ratioed_image.shape[0],ratioed_image.shape[1])))

        
    column_filtered_Fe_oxide_classification = []
    for i in range(len(reshaped_thresholded_a_530)):
        for j in range(len(reshaped_thresholded_a_530[i])):  
            if reshaped_thresholded_a_530[i][j]==1 and reshaped_thresholded_a_530[i][j+1]==0 and reshaped_thresholded_a_530[i][j-1]==0:
                column_filtered_Fe_oxide_classification.append(0)
            elif reshaped_a_530[i][j] < BD625_reshaped[i][j]:
                column_filtered_Fe_oxide_classification.append(0)
            else:
                column_filtered_Fe_oxide_classification.append(reshaped_thresholded_a_530[i][j])
    reshaped_column_filtered_Fe_oxide_classification = np.reshape(np.array(column_filtered_Fe_oxide_classification), (ratioed_image.shape[0],ratioed_image.shape[1]))

    path = #path to save data
    cv2.imwrite(path, reshaped_column_filtered_Fe_oxide_classification)
    output_folder = #path to save data
    base_name = str(image_id) + str('_column_filtered_Fe_oxide_classification') 
    header_file = os.path.join(output_folder, base_name + '.hdr')
    spectral.envi.save_image(header_file, reshaped_column_filtered_Fe_oxide_classification, interleave='bil', dytpe=np.float32, metadata = metadata_dict)
    path = #path to save data
    matplotlib.image.imsave(path, reshaped_column_filtered_Fe_oxide_classification, vmin=0, vmax=1, cmap='gray')
    
    
    column_filtered_thresholded_gauss_a = []
    for i in range(len(reshaped_thresholded_a_530)):
        for j in range(len(reshaped_thresholded_a_530[i])):
            if reshaped_column_filtered_Fe_oxide_classification[i][j]==1:
                column_filtered_thresholded_gauss_a.append(reshaped_a_530[i][j])
            else:
                column_filtered_thresholded_gauss_a.append(0)
    reshaped_column_filtered_thresholded_gauss_a = np.reshape(np.array(column_filtered_thresholded_gauss_a), (ratioed_image.shape[0],ratioed_image.shape[1]))
    path = #path to save data
    cv2.imwrite(path, reshaped_column_filtered_thresholded_gauss_a)
    output_folder = '#path to save data
    base_name = str(image_id) + str('_Fe_oxide_column_filtered_thresholded_gauss_a') 
    header_file = os.path.join(output_folder, base_name + '.hdr')
    spectral.envi.save_image(header_file, reshaped_column_filtered_thresholded_gauss_a, interleave='bil', dytpe=np.float32, metadata = metadata_dict)
    path = #path to save data
    matplotlib.image.imsave(path, reshaped_column_filtered_thresholded_gauss_a, vmin=0, vmax=0.03, cmap='gray')
    
    
    path = #path to save data
    np.savetxt(path, a_530)
    path = #path to save data
    np.savetxt(path, b_530)
    path = #path to save data
    np.savetxt(path, c_530)
    path = #path to save data
    np.savetxt(path, straight_line_cont_rem_530)


# In[ ]:




