% generate an ensemble of different kernel
% used to determine best kernel out of selection
function output = choose_from_kernel_ensemble(p, image)

    import engines.debluring.*
    
    mode_list = {'gaussian','gaussian_asym','diagonal','boxblur'};
    params_list = {};

    xvalues = 1 : 1 : 2;
    kernel_size = p.kernel_size;
    %gaussian kernel
    params_list{end+1} = {kernel_size,xvalues};
    %gaussian_asym kernel
    [X,Y] = meshgrid(xvalues,xvalues);
    params_list{end+1} = {kernel_size,[X(:),Y(:)].'};
    %diagonal kernel
    weights = 0 : 0.1 : 1;
    params_list{end+1} = {kernel_size,weights};
    %boxblur kernel
    params_list{end+1} = {kernel_size,[1]};
    
    settings = struct();
    settings.SNRt = p.SNRt;                     % SNRt = 0.2071 for 1/2 bit threshold for average of 2 images
                                                % SNRt = 0.5 for 1 bit threshold for average of 2 images
    settings.thickring = p.thickring;           % thick ring in Fourier domain
    settings.auto_binning = p.auto_binning;     % bin FRC before calculating rings, it makes calculations faster 
    settings.max_rings = p.max_rings;           % maximal number of rings if autobinning is used 
    settings.freq_thr = p.freq_thr;             % mimimal freq value where resolution is detected  
    settings.pixel_size = p.dx_spec;            % size of pixel in angstrom
    settings.mask = ones(size(image,1),size(image,2));
    settings.lucy_iters = p.deconvlucy_iters;
    settings.damping_threshold = p.damping_threshold;
    settings.correlation_threshold = p.correlation_threshold;
    
    max_res = zeros(1,4);
    max_idx = zeros(1,4);
    for mode_idx = 1 : 1 : 4
        p.kernel_type = mode_list{mode_idx};
        kernel_size = params_list{mode_idx}{1};
        param_set = params_list{mode_idx}{2};
        integrated_res = zeros(1,size(param_set,2));               
        for param_idx = 1 : 1 : size(params_list{mode_idx}{2},2)          
            kernel_fcn = @(params) init_kernel(p,kernel_size, params);
            integrated_res(1,param_idx) = FSC_3d_e(kernel_fcn(param_set(:,param_idx)),image,settings); % integrated maximal resolution
        end
        [max_res(mode_idx),max_idx(mode_idx)] = max(integrated_res,[],'all','linear');
    end
    [~,idx_highest_res] = max(max_res,[],'all','linear');

    params = params_list{idx_highest_res}{max_idx(idx_highest_res)}; 
    optimal_kernel = init_kernel(p,kernel_size, params);

    output = {optimal_kernel,params};    
end    