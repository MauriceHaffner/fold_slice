function mask = kernel_to_mask(target_size,kernel)
    
    % Basically the same as the Fourier Transform of kernel, i.e. projection to the farfield
    import engines.debluring.*

    padded_kernel = pad_kernel(target_size,kernel);
    mask = fftshift(fft2(padded_kernel));
        
end