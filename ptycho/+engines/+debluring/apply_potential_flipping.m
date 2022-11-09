% Apply potential flipping to the imaginary part of the object for pixel
% that fall below a threshold value. Aimed to increase separability of close
% features in the object.
function self = apply_potential_flipping(self, par, cache, iter)

    import engines.debluring.*
    import engines.GPU_MS.GPU_wrapper.*
    import math.*
    import utils.*
    import io.*

    object = Ggather(self.object{1});
    N_obj = size(object);
    object = zeros(N_obj(1),N_obj(2),par.Nlayers);
    laplace_kernel = [-1,-1,-1;-1,8,-1;-1,-1,-1];
    for ll = 1 : par.Nlayers
        object(:,:,ll) = Ggather(self.object{ll});
    end
    object_imag = imag(object);
    cell_of_GPU_arrays = {};
    for ll = 1 : par.Nlayers
        mean_imag = mean(object(cache.objectROI{:},ll),'all');
        flip_idx = find(object_imag(:,:,ll) < mean_imag * par.p.flip_threshold & object_imag(:,:,ll) > 0);
        object(flip_idx) = object(flip_idx,ll) - 2 * 1j * object_imag(flip_idx,ll);
        cell_of_GPU_arrays{end+1} = Garray(object(:,:,ll));
    end  
    self.object = cell_of_GPU_arrays;   

end