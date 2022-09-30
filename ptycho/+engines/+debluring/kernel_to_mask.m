function mask = kernel_to_mask(target_size,kernel)
    
    import engines.debluring.*

    padded_kernel = pad_kernel(target_size,kernel);
    mask = fft2(fftshift(padded_kernel));
        
end