function kernel = arbitrary_kernel(size,param)

kernel = reshape(param, [size,size]);

% kernel = kernel / sum(kernel,'all'); <--- specify this for contsrained optimization by Aeq = ones(lenght(params)), beq = 1

end