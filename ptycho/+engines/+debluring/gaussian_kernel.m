function kernel = gaussian_kernel(size,params)

    sigma = params(1);
    
    if (mod(size,2)~=0)
        size = size + 1;
    end
    x = -1/2:1/(size-1):1/2;
    [X,Y] = meshgrid(x,x);
    kernel = exp( -(X.^2+Y.^2)/(2*sigma^2) );
    kernel = kernel / sum(kernel, 'all');
    
end