function kernel = diagonal_kernel(size)
    
    if (mod(size,2)~=0)
        size = size + 1;
    end
    x = 0 : 1 : size-1;
    x = x - (size-1)/2;
    kernel = abs(x);
    kernel((n-1)/2 + 1) = size;
    kernel = 1 ./ kernel;
    kernel = kernel / sum(kernel,"all");
    kernel = diag(kernel);
        
end