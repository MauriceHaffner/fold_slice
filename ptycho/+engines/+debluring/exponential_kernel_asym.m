function kernel = exponential_kernel_asym(size,params)

    % Make a asymmetric in x and y exponentially decaying kernel
    
    sigma_x = params(1);
    sigma_y = params(2);
    a = params(3);

    if (mod(size,2)~=0)
        size = size + 1;
    end
    x = -1/2:1/(size-1):1/2;
    [X,Y] = meshgrid(x,x);
    kernel = exp( -abs(X).^a/(2*abs(sigma_x)) -abs(Y).^a/(2*abs(sigma_y)) );
    kernel = kernel / sum(kernel, 'all');
    
end