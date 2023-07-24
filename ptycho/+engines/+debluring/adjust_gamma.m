function self = adjust_gamma(self,par,cache)

    % Apply gamma adjustment to the object for each slicee

    import engines.debluring.*
    import engines.GPU_MS.GPU_wrapper.*
    import math.*
    import utils.*
    import io.*

    % obtain the object
    object = Ggather(self.object{1});
    object = zeros(size(object,1),size(object,2),par.Nlayers);
    cell_of_GPU_arrays = {};
    for ll = 1 : par.Nlayers
        object(:,:,ll) = Ggather(self.object{ll});
    end
           
   % Apply the adjustment 
    for idx = 1 : 1 : par.p.Nlayers
        Im = imadjust(mat2gray(angle(object(:,:,idx))),[],[],par.p.gamma);
        max_angle = max(angle(object(:,:,idx)),[],"all");
        min_angle = min(angle(object(:,:,idx)),[],"all");
        % Put the object back to GPU        
        cell_of_GPU_arrays{end+1} = Garray(abs(object(:,:,idx)) .* exp (1j * (Im*(max_angle-min_angle) + min_angle)));
    end  
    self.object = cell_of_GPU_arrays;
    
end