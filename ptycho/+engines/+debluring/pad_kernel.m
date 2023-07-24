function padded_kernel = pad_kernel(target_size,kernel)

    % Pad a kernel with zeros
    if (ndims(kernel) == 2)
        padding_size = floor((target_size - size(kernel))/2);
        padded_kernel = padarray(kernel, [padding_size(1), padding_size(2)], 0, 'both');
        padding_size = target_size-size(padded_kernel);
        padded_kernel = padarray(padded_kernel, [padding_size(1), padding_size(2)],0, 'pre');
    elseif (ndims(kernel) ==  3)
        padding_size = floor((target_size - [size(kernel,1),size(kernel,2)])/2);
        padded_kernel = padarray(kernel, [padding_size(1), padding_size(2),0], 0, 'both');
        padding_size = target_size-[size(padded_kernel,1),size(padded_kernel,2)];
        padded_kernel = padarray(padded_kernel, [padding_size(1), padding_size(2),0],0, 'pre');
    else
        error('Invalid kernel size.')
    end           
end