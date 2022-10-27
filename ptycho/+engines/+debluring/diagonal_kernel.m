function kernel = diagonal_kernel(size,params)
   
    import engines.debluring.*
    weight= params(1);
    if (mod(size,2)==0)
        size = size + 1;
    end
    kernel = weight * horizontal_kernel(size) + (1-weight) * horizontal_kernel(size).';
    kernel = kernel / sum(kernel,'all');
        
end