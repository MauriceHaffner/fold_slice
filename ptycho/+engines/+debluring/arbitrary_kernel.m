function kernel = arbitrary_kernel(size,param)

    % Reshape a an arbitray paramter vector to the desired kernel size
    kernel = reshape(param, [size,size]);

end