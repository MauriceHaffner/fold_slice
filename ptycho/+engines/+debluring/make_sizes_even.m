function cropped_image = make_sizes_even(image)

    % Crop a image so that all dimension sizes are even
    
    original_size = size(image);    
    target_size = [2*floor(original_size(1)/2),2*floor(original_size(2)/2)];
    cropped_image = image(1:target_size(1),1:target_size(2),:);
   
end