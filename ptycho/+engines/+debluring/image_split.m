function output = image_split(image)

    image_size = size(image);
    sub_image_size = image_size;
    sub_image_size(2:end) = ceil(sub_image_size(2:end)/2);
    
    img1 = zeros(sub_image_size);
    img2 = zeros(sub_image_size);
    img3 = zeros(sub_image_size);
    img4 = zeros(sub_image_size);

    img1 = image(1:2:end,1:2:end,:); % odd-odd split
    img2 = image(2:2:end,2:2:end,:); % even-even split
    img3 = image(1:2:end,2:2:end,:); % odd-even split
    img4 = image(2:2:end,1:2:end,:); % even-odd-split

    output = {img1,img2,img3,img4}; % return multiple multidimensional arrays as cell array

end