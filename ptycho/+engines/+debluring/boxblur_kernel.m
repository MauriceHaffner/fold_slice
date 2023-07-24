function kernel = boxblur_kernel(size)

    % Kernel with every pixel having same weight
    kernel = 1/size^2 * ones(size,size);

end