function kernel = gaussian_kernel_asym(size,params)

    sigma_x = params(1);
    sigma_y = params(2);

    if (mod(size,2)~=0)
        size = size + 1;
    end
    x = -1/2:1/(size-1):1/2;
    [X,Y] = meshgrid(x,x);
    kernel = exp( -X.^2/(2*sigma_x^2) -Y.^2/(2*sigma_y^2));
    kernel = kernel / sum(kernel, 'all');
    
end