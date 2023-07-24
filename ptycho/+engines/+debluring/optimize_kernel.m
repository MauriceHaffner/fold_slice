
function output = optimize_kernel(p, image)

    import engines.debluring.*
    
    if (strcmp(p.kernel_type,'horizontal') || strcmp(p.kernel_type,'vertical') || strcmp(p.kernel_type,'diagonal'))
        warning('Specified kernel type is not optimizable.')
        optimal_kernel = kernel_fcn(x0);
        return
    end

    % settings = struct();
    % settings.SNRt = p.SNRt;                     % SNRt = 0.2071 for 1/2 bit threshold for average of 2 images
    %                                             % SNRt = 0.5 for 1 bit threshold for average of 2 images
    % settings.thickring = p.thickring;           % thick ring in Fourier domain
    % settings.auto_binning = p.auto_binning;     % bin FRC before calculating rings, it makes calculations faster 
    % settings.max_rings = p.max_rings;           % maximal number of rings if autobinning is used 
    % settings.freq_thr = p.freq_thr;             % mimimal freq value where resolution is detected  
    % settings.pixel_size = p.dx_spec;            % size of pixel in angstrom
    % settings.mask = ones([size(image,1),size(image,2)]);
    % settings.lucy_iters = p.deconvlucy_iters;
    % settings.damping_threshold = p.damping_threshold;
    % settings.correlation_threshold = p.correlation_threshold;
    
    check = false;
    storage_dim = [p.kernel_size,p.kernel_size,size(p.kernel_params,1)];

    kernel_fcn = @(params) init_kernel(p,p.kernel_size, params);    
    %optimal_kernels = zeros(storage_dim);
    %x = zeros(size(p.kernel_params,1),size(p.kernel_params,2));
    loadvars = {"p","object"};
    %kernel_L2_errors = zeros(1,size(p.kernel_params,1));
    initial_sparsity = sum(log10(abs(fft2(imag(image))) ./ max(abs(fft2(imag(image))),[],"all"))> -1,"all")
    for idx=1:1:size(p.kernel_params,1)
        x0 = p.kernel_params(idx,:);        

        [image_deconv,psfr] = deconvblind(mat2gray(imag(image)),kernel_fcn(x0)+rand(p.kernel_size,p.kernel_size)*1e-2);
        [x_opt,L2_error] = minimize_kernel_diff_L2(psfr,kernel_fcn,x0);
        sparsity = sum(log10(abs(fft2(image_deconv)) ./ max(abs(fft2(image_deconv)),[],"all"))> -1,"all");
        %keyboard
        %x(idx,:) = x_opt;
        %kernel_L2_errors(1,idx) = L2_error;
        %optimal_kernels(:,:,idx) = psfr;
        % if (L2_error < p.target_error) break;
        % end
        if (sparsity < p.target_sparsity * initial_sparsity) break;
        end
    end
    %[least_error,optimal_idx] = min(kernel_L2_errors);
    sparsity
    %output = {optimal_kernels(:,:,optimal_idx),x(optimal_idx,:),least_error};
    output = {psfr,x_opt,L2_error};  
end