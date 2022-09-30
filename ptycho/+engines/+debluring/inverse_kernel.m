function kernel_inv = inverse_kernel(p,target_size,kernel)
    
    import engines.debluring.*

    padded_kernel = pad_kernel(target_size,kernel);
    fft_kernel = fft2(fftshift(padded_kernel));
    kernel_inv = conj(fft_kernel) ./ abs(fft_kernel);
    kernel_inv(isnan(kernel_inv)) = 0;
    kernel_inv(isinf(kernel_inv)) = 0;
    if isfield(p, 'smooth_kernel') && p.smooth_kernel && isfield(p, 'smooth_kernel_width')
    gauss_kernel = gaussian_kernel(9,p.smooth_kernel_width);
    kernel_inv = conv2(kernel_inv,gauss_kernel,"same");
    end
        
end