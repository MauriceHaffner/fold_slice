function kernel = gaussian_kernel(size,params)

    % Make a asymmetric in x and y gaussian kernel
    
    sigma = params(1);
    
    if (mod(size,2)~=1)
        size = size + 1;
    end
    %x = -1/2/2:1/(size-1):1/2;
    x = -(size-1)/2 : 1 : (size-1)/2;
    [X,Y] = meshgrid(x,x);
    kernel = exp( -(X.^2+Y.^2)/(2*sigma^2) );
    kernel = kernel / sum(kernel, 'all');
    
end