function padded_kernel = pad_kernel(target_size,kernel)

    %if mod(target_size(1),2) | mod(target_size(2),2) % one or both dimensions are uneven
    padding_size = floor((target_size - size(kernel))/2);
    padded_kernel = padarray(kernel, [padding_size(1), padding_size(2)], 0, 'both');
    padding_size = target_size-size(padded_kernel);%[mod(target_size(1),2),mod(target_size(2),2)];
    padded_kernel = padarray(padded_kernel, [padding_size(1), padding_size(2)], 'pre');
    % else   % both dimensions are even 
    %     padding_size = (target_size - size(kernel))/2;
    %     keyboard
    %     padded_kernel = padarray(kernel, [padding_size(1), padding_size(2)], 0, 'both');
    % end

end