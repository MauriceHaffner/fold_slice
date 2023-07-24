% DEPRECATED
% function self = apply_sharpening(self, par, cache, iter)

%     % Apply sharpening to the object by increasing the gradient of the imaginary part
%     import engines.debluring.*
%     import engines.GPU_MS.GPU_wrapper.*
%     import math.*
%     import utils.*
%     import io.*

%     object = Ggather(self.object{1});
%     N_obj = size(object);
%     object = zeros(N_obj(1),N_obj(2),par.Nlayers);
%     laplace_kernel = [0,-1,0;-1,4,-1;-0,-1,0];
%     for ll = 1 : par.Nlayers
%         object(:,:,ll) = Ggather(self.object{ll});
%     end
%     object_real = real(object);
%     object_imag = imag(object);
%     cell_of_GPU_arrays = {};
%     % Compute the second derivative of the imaginary part of the object
%     % and add it to the object. Probably working only when WPOA holds.
%     for ll = 1 : par.Nlayers
%         imag_der = conv2(conv2(object_imag(:,:,ll),laplace_kernel,'same'),laplace_kernel,'same');
%         cell_of_GPU_arrays{end+1} = Garray(object(:,:,ll) + par.p.sharpening_weight * 1j * imag_der);
%     end  
%     self.object = cell_of_GPU_arrays;   

% end    