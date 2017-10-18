function [ extrapolated_image ] = crimminisiExtrapolate( image_in, image_bounds, patch_size )
    %Create the mask
    mask = ones(size(image_in,1), size(image_in,2));
    mask(image_bounds(1,1):image_bounds(1,2),...
        image_bounds(2,1):image_bounds(2,2)) = 0;
    mask = logical(mask);
    
    %Fill the missing parts
    [extrapolated_image,~,~,~] = inpainting(image_in,mask,patch_size);
end

