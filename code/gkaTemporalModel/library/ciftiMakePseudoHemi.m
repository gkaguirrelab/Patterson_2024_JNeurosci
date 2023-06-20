function ciftiOut = ciftiMakePseudoHemi(ciftiIn)

    % Copy input to output
    ciftiOut = ciftiIn;

    % Get the left and right hemisphere data
    leftHemi = cifti_struct_dense_extract_surface_data(ciftiIn, 'CORTEX_LEFT');
    rightHemi = cifti_struct_dense_extract_surface_data(ciftiIn, 'CORTEX_RIGHT');
    
    % Average
    averageHemi = nanmean([rightHemi,leftHemi], 2);
    
    % Replace the left hemi and right hemi values with the new average
    ciftiOut = cifti_struct_dense_replace_surface_data(ciftiOut, averageHemi, 'CORTEX_LEFT');
    ciftiOut = cifti_struct_dense_replace_surface_data(ciftiOut, averageHemi, 'CORTEX_RIGHT');

end