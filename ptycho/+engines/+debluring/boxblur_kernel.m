function kernel = boxblur_kernel(size)

    kernel = 1/size^2 * ones(size,size);

end