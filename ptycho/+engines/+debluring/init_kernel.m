% call custom functions for debluring/ deconvolution

function output_kernel = init_kernel(p, params, varargin)

    import engines.debluring.*

    size = p.kernel_size;
    switch p.kernel_type

    case 'gaussian'
        kernel = gaussian_kernel(size,params);
    case 'gaussian_asym'
        kernel = gaussian_kernel_asym(size,params);
    case 'arbitrary'
        kernel = arbitrary_kernel(size,params);
    case 'horizontal'
        kernel = horizontal_kernel(size,params);
    case 'vertical'
        kernel = transpose(horizontal_kernel(size,params));
    case 'diagonal'
        kernel = diagonal_kernel(size);
    otherwise
        try
            custom_fcn = str2func(type);
            kernel = custom_fcn(size,params);
        catch
            error('Unspecified kernel type or not found.');
        end
    end
        
    if isfield(p, 'kernel_residual_style') && p.kernel_residual_style
    % allow to optimize for a kernel, that is roughly equal to the identity convolution kernel    
        identity_kernel = zeros(size);
        identity_kernel( ceil(size/2), ceil(size/2))  = 1;
        kernel = identity_kernel + params(end) * kernel;
        kernel = kernel / sum(kernel, 'all');
    end

    if strcmp(p.output_kernel, 'inv')
        target_size = varargin;
        output_kernel = inverse_kernel(p,target_size,kernel);
    elseif strcmp(p.output_kernel, 'real')
        output_kernel = kernel;
    end        

end