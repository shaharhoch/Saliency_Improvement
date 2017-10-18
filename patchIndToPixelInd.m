function [ pixel_index ] = patchIndToPixelInd( patch_index, patch_size, image )
% Returns the index of the top left corner of the patch
    num_patches_col = floor(size(image,1)/patch_size);
    patch_index_col = mod(patch_index-1, num_patches_col)+1;
    patch_index_row = floor((patch_index-1)/num_patches_col)+1;
    
    pixel_index = ([patch_index_col, patch_index_row]-1)*patch_size+1;

end

