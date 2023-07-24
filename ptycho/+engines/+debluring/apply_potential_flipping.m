function [self,par,cache] = apply_potential_flipping(self, par, cache, iter)

    % Apply potential flipping to the imaginary part of the object for pixel
    % that fall below a threshold value. Aimed to increase separability of close
    % features in the object.

    import engines.debluring.*
    import engines.GPU_MS.GPU_wrapper.*
    import math.*
    import utils.*
    import io.*

    % Obtain the object from the self struct
    object = Ggather(self.object{1});
    N_obj = size(object);
    object = zeros(N_obj(1),N_obj(2),par.Nlayers);
    for ll = 1 : par.Nlayers
        object(:,:,ll) = Ggather(self.object{ll});
    end    
    cell_of_GPU_arrays = {};
    for ll = 1 : par.Nlayers
        object_imag = imag(object(:,:,ll));
        max_imag = max(object_imag,[],'all');
        flip_idx = find((object_imag < 0.3 * max_imag)); % Find the indices where to flip. 30 % of maximum imaginary part was found empirically to work well.
        object_imag(flip_idx) = -1 * object_imag(flip_idx);
        cell_of_GPU_arrays{end+1} = Garray(real(object(:,:,ll)) + 1j * object_imag);
    end  
    self.object = cell_of_GPU_arrays;   

end