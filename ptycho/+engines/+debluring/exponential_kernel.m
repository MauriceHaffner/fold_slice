function kernel = exponential_kernel(size,params)

    % Make a symmetric in x and y exponentially decaying kernel
    
    sigma = params(1);
    a = params(2);

    if (mod(size,2)~=0)
        size = size + 1;
    end
    x = -1/2:1/(size-1):1/2;
    [X,Y] = meshgrid(x,x);
    kernel = exp( -(abs(X).^a + abs(Y).^a)/(2*abs(sigma)));
    kernel = kernel / sum(kernel, 'all');
    
end