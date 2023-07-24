function kernel = horizontal_kernel(size)

    % Make a kernel with only horizontal non-zero elements
    
    if (mod(size,2)==0)
        size = size + 1;
    end
    kernel = zeros(size,size);
    kernel((size-1)/2,:) = 1/size * ones(1,size);
        
end