function [ closest_patch_ind ] = getClosesPatchInd( patch_ind, possible_patches, image, patch_size )
    min_dist = intmax;
    closest_patch_ind = 0;
    for i=1:length(possible_patches)
        dist = getPatchDistance( patch_ind, possible_patches(i), image, patch_size );
        if(dist < min_dist)
            min_dist = dist; 
            closest_patch_ind = possible_patches(i);
        end
    end
    assert(closest_patch_ind ~= 0);
end

function [ dist ] = getPatchDistance( patch1_ind, patch2_ind, image, patch_size )
    pixel_index1 = patchIndToPixelInd( patch1_ind, patch_size, image );
    pixel_index2 = patchIndToPixelInd( patch2_ind, patch_size, image );
    
    patch1 = image(pixel_index1(1):pixel_index1(1)+patch_size-1,...
        pixel_index1(2):pixel_index1(2)+patch_size-1,:);
    patch2 = image(pixel_index2(1):pixel_index2(1)+patch_size-1,...
        pixel_index2(2):pixel_index2(2)+patch_size-1,:);
    
    diff = patch1(:)-patch2(:);
    dist = norm(diff);
end
