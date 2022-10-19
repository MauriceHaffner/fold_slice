function self = apply_deconvolution(self, par, cache, iter)

    import engines.debluring.*
    import engines.GPU_MS.GPU_wrapper.*
    import math.*
    import utils.*
    import io.*

    % obtain the object ROI and apply deconvolution and obtain the optimal kernel
    objectROI = Ggather(self.object{1}(cache.object_ROI{:},:));
    N_obj_roi = size(objectROI);
    objectROI = zeros(N_obj_roi(1),N_obj_roi(2),par.Nlayers);
    for ll = 1 : par.Nlayers
        objectROI(:,:,ll) = Ggather(self.object{ll}(cache.object_ROI{:},:));
    end
    objectROI = make_sizes_even(objectROI);
    output_cell = optimize_kernel(par.p,objectROI);
    par.p.optimal_kernel = output_cell{1};
    par.p.optimal_kernel_params = output_cell{2};    
    % now apply the deconvolution to the whole object
    %if par.p.deconv_final_iter == iter
    object = Ggather(self.object{1});
    N_obj = size(object);
    object = zeros(N_obj(1),N_obj(2),par.Nlayers);
    for ll = 1 : par.Nlayers
    object(:,:,ll) = Ggather(self.object{ll});
    end
    object_real = real(object);
    object_imag = imag(object);
    object_imag_deconv = deconvlucy(object_imag, par.p.optimal_kernel,par.p.deconvlucy_iters,par.p.SNRt);
    cell_of_GPU_arrays = {};
    for ll = 1 : par.Nlayers
        cell_of_GPU_arrays{end+1} = Garray((1-par.p.deconv_weight)*(object_real(:,:,ll) + 1j * object_imag(:,:,ll)) + par.p.deconv_weight*(object_real(:,:,ll) + 1j * object_imag_deconv(:,:,ll)));
    end  
    self.object = cell_of_GPU_arrays;
        %return; % skip the update of the mask if arrived at the final iteration for deconvolution
    %end    
    % update the mask, used within next engine to counter effect of convolution noise
    mask = par.p.fmask;
    mask_size = [size(mask,1),size(mask,2)];
    kernel_mask = kernel_to_mask(mask_size,par.p.optimal_kernel);
    binary_mask = abs(kernel_mask) > par.p.mask_threshold * max(abs(kernel_mask),[],'all');
    par.p.fmask = logical(binary_mask); %.* mask

    % mask = par.p.fmask;
    % mask_size = [size(mask,1),size(mask,2)];
    % par.p.fmask = kernel_to_mask(mask_size,par.p.optimal_kernel) .* mask;
    
end    

