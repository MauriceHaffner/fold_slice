
function output = optimize_kernel(p, image)

    import engines.debluring.*
    
    if (strcmp(p.kernel_type,'horizontal') || strcmp(p.kernel_type,'vertical') || strcmp(p.kernel_type,'diagonal'))
        warning('Specified kernel type is not optimizable. \n')
        optimal_kernel = kernel_fcn(x0);
        return
    end

    settings = struct();
    settings.SNRt = p.SNRt;                     % SNRt = 0.2071 for 1/2 bit threshold for average of 2 images
                                                % SNRt = 0.5 for 1 bit threshold for average of 2 images
    settings.thickring = p.thickring;           % thick ring in Fourier domain
    settings.auto_binning = p.auto_binning;     % bin FRC before calculating rings, it makes calculations faster 
    settings.max_rings = p.max_rings;           % maximal number of rings if autobinning is used 
    settings.freq_thr = p.freq_thr;             % mimimal freq value where resolution is detected  
    settings.pixel_size = p.dx_spec;            % size of pixel in angstrom
    settings.mask = ones([size(image,1),size(image,2)]);
    settings.lucy_iters = p.deconvlucy_iters;
    settings.damping_threshold = p.damping_threshold;
    settings.correlation_threshold = p.correlation_threshold;
    
    iter_max = 1;
    check = false;
    storage_dim = [p.kernel_size,p.kernel_size,1];
    if isfield(p, 'multislice_deconvolution') && (p.Nlayers) > 1 && p.multislice_deconvolution
        iter_max = p.Nlayers;
        storage_dim = [p,kernel_size,p.kernel_size,p.Nlayers];
        check = true;   
    end

    kernel_fcn = @(params) init_kernel(p, params);    
    options = optimset('MaxFunEvals',p.FSC_evals,'MaxIter', p.fminsearch_evals, 'TolFun', p.resolution_tol);
    optimal_kernel = zeros(storage_dim);
    x = zeros(iter_max,size(p.kernel_params,2));
    fvals = zeros(iter_max,1);
    exitflags = zeros(iter_max,1); 
    for idx=1:1:iter_max
        x0 = p.kernel_params(idx,:);
        if check
            opt_fcn = @(x) (-1) * FSC_3d_e(kernel_fcn(x),image(:,:,idx),settings); % = (-1) * maximal resolution
        else
            opt_fcn = @(x) (-1) * FSC_3d_e(kernel_fcn(x),image,settings); % = (-1) * maximal resolution
        end
        A =         p. A_linear_ineq_constr;
        b =         p. b_linear_ineq_contsr;
        Aeq =       p. A_linear_eq_constr;
        beq =       p. b_linear_eq_contsr;
        lb =        p. kernel_params_lb;
        ub =        p. kernel_params_ub;
        nonclon =   p. nonlinear_constr;

        try
            [x(idx,:), fvals(idx), exitflags(idx)] = fmincon(opt_fcn,x0,A,b,Aeq,beq,lb,ub,nonclon,options);
            optimal_kernel(:,:,idx)= kernel_fcn(x);
        catch
            warning('The deconvolution algorithm failed due to an error with fmincon. Set optimal kernel to identity kernel. \n');
            identity_kernel = zeros(storage_dim);
            identity_kernel(ceil(storage_dim(1)/2),ceil(storage_dim(2)/2)) = 1;
            optimal_kernel(:,:,idx) = identity_kernel;
            x = [1];
        end    
    end

    if ~all(exitflags)
        warning('The deconvolution algorithm did not converge at least once. Maybe check the settings. \n');
        keyboard
    end  

    output = {optimal_kernel,x};
      
end