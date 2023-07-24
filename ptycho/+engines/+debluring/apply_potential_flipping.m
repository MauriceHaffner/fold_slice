% Apply potential flipping to the imaginary part of the object for pixel
% that fall below a threshold value. Aimed to increase separability of close
% features in the object.
function [self,par,cache] = apply_potential_flipping(self, par, cache, iter)

    import engines.debluring.*
    import engines.GPU_MS.GPU_wrapper.*
    import math.*
    import utils.*
    import io.*

    object = Ggather(self.object{1});
    N_obj = size(object);
    object = zeros(N_obj(1),N_obj(2),par.Nlayers);
    for ll = 1 : par.Nlayers
        object(:,:,ll) = Ggather(self.object{ll});
    end    
    cell_of_GPU_arrays = {};
    for ll = 1 : par.Nlayers
        object_imag = imag(object(:,:,ll));
        std_imag = std(object_imag,0,"all");
        mean_imag = mean(object_imag,"all");
        max_imag = max(object_imag,[],'all');
        flip_idx = find((object_imag < 0.3 * max_imag));%max_imag * par.p.flip_threshold|| (object_imag > mean_imag + std_imag)- 0.8 * std_imag
        object_imag(flip_idx) = -1 * object_imag(flip_idx);
        cell_of_GPU_arrays{end+1} = Garray(object(:,:,ll) + 1j * object_imag);
    end  
    self.object = cell_of_GPU_arrays;   

end