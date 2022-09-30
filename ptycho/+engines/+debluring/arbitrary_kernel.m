function kernel = arbitrary_kernel(size,param)

kernel = reshape(param, [size,size]);
%kernel = kernel / sum(kernel,'all');

end