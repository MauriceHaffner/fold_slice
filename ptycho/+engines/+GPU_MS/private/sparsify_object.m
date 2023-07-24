% sparsify the object in the vertical (z) direction
function self = sparsify_object(self,par,cache)

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
depth = par.p.sparsify_depth;
% calculate a mask that is one only for values that are the maximum within a certain depth
sparsity_mask = ones(N_obj(1),N_obj(2),par.Nlayers);
for ll = 1 : depth : (par.Nlayers-depth+1)
    sparsity_mask(:,:,ll:ll+depth-1) = bsxfun(@ge,angle(object(:,:,ll:ll+depth-1)),max(angle(object(:,:,ll:ll+depth-1)),[],3));    
end
% sparsify the object by setting the phase to zero of those pixels who dont have the maximum phase withion the specified range
object_sparse = real(object) + 1j * imag(object) .* sparsity_mask;
% optionaly calculate the weighted average for smoother correction
if (par.p.sparsify_weight ~= 1) 
    object_sparse = (1-par.p.sparsify_weight) * object + par.p.sparsify_weight * object_sparse;
end

cell_of_GPU_arrays = {};
for ll = 1 : par.Nlayers
    cell_of_GPU_arrays{end+1} = Garray(object_sparse(:,:,ll));
end  
self.object = cell_of_GPU_arrays;

end