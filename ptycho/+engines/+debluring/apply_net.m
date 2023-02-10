function self = apply_net(self,par,cache)

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
        
    Im = zeros(size(object,1),size(object,2),3);
    residual = zeros(size(object,1),size(object,2),1);
    for idx = 1 : 1 : par.Nlayers

        Im(:,:,1) = mat2gray(angle(object(:,:,idx)));
        Im(:,:,2) = imadjust(Im(:,:,1),[],[],par.p.gamma);
        Im(:,:,3) = imbinarize(Im(:,:,1));
        % calculate the residual for all tiles
        for idx1 = 1 : par.p.patchSize : size(object,1)-par.p.patchSize
            for idx2 = 1 : par.p.patchSize : size(object,2)-par.p.patchSize
                residual(idx1:idx1+par.p.patchSize-1,idx2:idx2+par.p.patchSize-1,1) = predict(par.p.net,Im(idx1:idx1+par.p.patchSize-1,idx2:idx2+par.p.patchSize-1,:));
            end
        end

        % special treatment for those tiles at the very end, because they may exceed the boundary of the object
        % however, hence these regions are likely to be outside the actual ROI the effect should be rather small
        if (mod(size(object,1),par.p.patchSize) ~= 0)     
            for idx1 = 1 : par.p.patchSize : size(object,1)-par.p.patchSize
                residual(idx1:idx1+par.p.patchSize-1,size(object,2)-par.p.patchSize+1:end,1) = predict(par.p.net,Im(idx1:idx1+par.p.patchSize-1,size(object,1)-par.p.patchSize+1:end,:));
            end
        end

        if (mod(size(object,2),par.p.patchSize) ~= 0)
            for idx2 = 1 : par.p.patchSize : size(object,2)-par.p.patchSize
                residual(size(object,1)-par.p.patchSize+1:end,idx2:idx2+par.p.patchSize-1,1) = predict(par.p.net,Im(size(object,1)-par.p.patchSize+1:end,idx2:idx2+par.p.patchSize-1,:));
            end
        end

        residual(size(object,1)-par.p.patchSize+1:end,size(object,2)-par.p.patchSize+1:end,1) = predict(par.p.net,Im(size(object,1)-par.p.patchSize+1:end,size(object,2)-par.p.patchSize+1:end,:));

        max_angle = max(angle(object(:,:,idx)),[],"all");
        min_angle = min(angle(object(:,:,idx)),[],"all");        
        cell_of_GPU_arrays{end+1} = Garray(abs(object(:,:,idx)) .* exp (1j * ((Im(:,:,1) + par.p.residual_weight * residual)*max_angle + min_angle)));  
        
    end
    self.object = cell_of_GPU_arrays;

end