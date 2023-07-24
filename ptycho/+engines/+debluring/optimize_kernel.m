
function output = optimize_kernel(p, image)

    % Optimize a given kernel.
    import engines.debluring.*
    
    if (strcmp(p.kernel_type,'horizontal') || strcmp(p.kernel_type,'vertical') || strcmp(p.kernel_type,'diagonal'))
        warning('Specified kernel type is not optimizable.')
        optimal_kernel = kernel_fcn(x0);
        return
    end
    
    check = false;
    storage_dim = [p.kernel_size,p.kernel_size,size(p.kernel_params,1)];
    kernel_fcn = @(params) init_kernel(p,p.kernel_size, params);    
    loadvars = {"p","object"};
    % Use the sparsity in the Fourier Plane of the object as measure of improvement
    initial_sparsity = sum(log10(abs(fft2(imag(image))) ./ max(abs(fft2(imag(image))),[],"all"))> -1,"all")
    for idx=1:1:size(p.kernel_params,1)
        x0 = p.kernel_params(idx,:);        

        [image_deconv,psfr] = deconvblind(mat2gray(imag(image)),kernel_fcn(x0)+rand(p.kernel_size,p.kernel_size)*1e-2);
        [x_opt,L2_error] = minimize_kernel_diff_L2(psfr,kernel_fcn,x0);
        sparsity = sum(log10(abs(fft2(image_deconv)) ./ max(abs(fft2(image_deconv)),[],"all"))> -1,"all");
        if (sparsity < p.target_sparsity * initial_sparsity) break; % Once the target sparsity is reached stop the optimization
        end
    end
    output = {psfr,x_opt,L2_error};

end