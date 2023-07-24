function filtered_image = apply_inv_kernel(p,image,kernel,kernel_space_type)

    % Apply the inverse of a given convolution kernel yo an image
    
    import engines.debluring.*
    
    iter_max = 1;
    check = false;
    % Check if the object is multislice and should become deconvolved slice-wise
    if isfield(p, 'multislice_deconvolution') && p.Nlayers > 1 && p.multislice_deconvolution
        iter_max = p.Nlayers;  
    end
    
    filtered_image = zeros(size(image));
    fft_image = fft2(image);
    % The given kernel may be in real space or Fourier space
    if strcmp(kernel_space_type,'real')
        for idx=1:1:iter_max
            inv_kernel = inverse_kernel(p,kernel(:,:,idx));
            if check
            % apply the inverse kernel to the corresponding slice    
                filtered_image(:,:,idx) = ifft2(inv_kernel .* fft_image(:,:,idx));
            else
            % apply the inverse kernel to all slices       
                filtered_image = ifft2(inv_kernel .* fft_image);
            end    
        end
    elseif strcmp(kernel_space_type,'fourier')        
        filtered_image = ifft2(kernel .* fft_image);
    else
        warning('Kernel space for deconvolution was undefined! Either real or fourier. Deconvolution was NOT applied. \n ');
        filtered_image = image;
    end                
end