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
    if strcmp(par.p.kernel_selection,'optimization')
        output_cell = optimize_kernel(par.p,objectROI);
    elseif strcmp(par.p.kernel_selection,'choose_best') 
        output_cell = choose_from_kernel_ensemble(par.p,objectROI);
    else
        warning('Invalid kernel selection mode.');
        keyboard    
    end  
    par.p.optimal_kernel = output_cell{1};
    par.p.optimal_kernel_params = output_cell{2};    
    % now apply the deconvolution to the whole object
    if par.p.deconvolve_object
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
            cell_of_GPU_arrays{end+1} = Garray((1-par.p.deconv_weight_object)*(object_real(:,:,ll) + 1j * object_imag(:,:,ll)) + par.p.deconv_weight_object*(object_real(:,:,ll) + 1j * object_imag_deconv(:,:,ll)));
        end  
        self.object = cell_of_GPU_arrays;
    end  

    % update the mask, used within next engine to counter effect of convolution noise
    mask = par.p.fmask;
    mask_size = [size(mask,1),size(mask,2)];
    kernel_mask = kernel_to_mask(mask_size,par.p.optimal_kernel);
    binary_mask = abs(kernel_mask) > par.p.mask_threshold * max(abs(kernel_mask),[],'all');
    par.p.fmask = logical(binary_mask);

    %probe correction
    if par.p.deconvolve_probe        
        par.p.optimal_kernel = output_cell{1};
        par.p.optimal_kernel_params = output_cell{2};    
        probe = Ggather(self.probe{1});
        probe_size = size(probe);
        probe = zeros(probe_size(1),probe_size(2),par.Nmodes,probe_size(4));
        for n = 1 : 1 : par.Nmodes
            probe(:,:,n,:) = Ggather(self.probe{n});
        end
        if strcmp(par.p.kernel_selection,'optimization')
            output_cell = optimize_kernel(par.p,probe(:,:,:,2));
        elseif strcmp(par.p.kernel_selection,'choose_best') 
            output_cell = choose_from_kernel_ensemble(par.p,probe(:,:,:,2));
        else
            warning('Invalid kernel selection mode.');
            keyboard    
        end
        par.p.optimal_kernel = output_cell{1};
        par.p.optimal_kernel_params = output_cell{2};      
        phase_real = real(probe(:,:,:,2));
        phase_imag = imag(probe(:,:,:,2));
        phase_imag_deconv = deconvlucy(phase_imag, par.p.optimal_kernel,par.p.deconvlucy_iters,par.p.SNRt);
        cell_of_GPU_arrays = {};
        probe_new = zeros(probe_size);
        for jj = 1 : 1 : par.Nmodes
            probe_new(:,:,1,1) = probe(:,:,jj,1);
            probe_new(:,:,1,2) = (1-par.p.deconv_weight_probe)*(phase_real(:,:,jj) + 1j * phase_imag(:,:,jj)) + par.p.deconv_weight_probe*(phase_real(:,:,jj) + 1j * phase_imag_deconv(:,:,jj));
            cell_of_GPU_arrays{end+1} = Garray(probe_new);          
        end  
        self.probe = cell_of_GPU_arrays;
    end
    
end    

