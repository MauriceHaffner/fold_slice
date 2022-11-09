function FFTkernel = fft_kernel(target_size,kernel)
    
    import engines.debluring.*

    padded_kernel = pad_kernel(target_size,kernel);
    FFTkernel = fft2(fftshift(padded_kernel));
        
end