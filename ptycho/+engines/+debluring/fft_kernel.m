function FFTkernel = fft_kernel(target_size,kernel)
    
    % Compute the Fourier Transform of a given kernel with padding added
    import engines.debluring.*

    padded_kernel = pad_kernel(target_size,kernel);
    FFTkernel = fft2(fftshift(padded_kernel));
        
end