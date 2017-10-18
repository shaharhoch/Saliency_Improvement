function [ replaced_patches_idx, left_patches_idx ] = getPatchesIdx( image, image_mask, patch_size, replace_ratio )
    num_of_patches = floor(size(image,1)/patch_size) * floor(size(image,2)/patch_size);
    
    ind_array = [];
    for i=1:num_of_patches
        if(~isPatchInsideMask( i, image_mask, patch_size ))
            ind_array = [ind_array, i];
        end
    end
    num_replaced = floor(replace_ratio*length(ind_array));
            
    replaced_patches_idx = datasample(ind_array, num_replaced, 'Replace', false);
    left_patches_idx = setdiff(ind_array, replaced_patches_idx);
end

function [ is_inside ] = isPatchInsideMask( patch_ind, mask, patch_size )
    is_inside = false;
    pixel_ind = patchIndToPixelInd(patch_ind, patch_size, mask);
    
    if(mask(pixel_ind(1), pixel_ind(2)) == 1)
        is_inside = true;
        return
    end
    
    if(mask(pixel_ind(1)+patch_size-1, pixel_ind(2)) == 1)
        is_inside = true;
        return
    end
    
    if(mask(pixel_ind(1), pixel_ind(2)+patch_size-1) == 1)
        is_inside = true;
        return
    end
    
    if(mask(pixel_ind(1)+patch_size-1, pixel_ind(2)+patch_size-1) == 1)
        is_inside = true;
        return
    end
end
    
    
    
    
    
    
