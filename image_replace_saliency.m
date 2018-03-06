function [ out_image ] = image_replace_saliency( in_img, object_mask, decimination_ratio, patch_size )
%IMAGE_REPLACE_SALIENCY Summary of this function goes here
%   Detailed explanation goes here

in_img_rgb = rgb2lab(in_img);
out_image_rgb = in_img_rgb; 

%Get replaced and left patches index
[replaced_patches_idx, left_patches_idx] = getPatchesIdx(in_img_rgb,...
    object_mask, patch_size, decimination_ratio );

% Replace patches
for ind=1:length(replaced_patches_idx)
    closest_patch_ind = getClosesPatchInd(replaced_patches_idx(ind),...
        left_patches_idx, in_img_rgb, patch_size);

    %Replace the patch
    orig_pixel_index = patchIndToPixelInd(replaced_patches_idx(ind), patch_size, out_image_rgb);
    replace_pixel_ind = patchIndToPixelInd(closest_patch_ind, patch_size, out_image_rgb);

    cur_patch = in_img_rgb(orig_pixel_index(1):orig_pixel_index(1)+patch_size-1,...
        orig_pixel_index(2):orig_pixel_index(2)+patch_size-1,:);
    replace_patch = in_img_rgb(replace_pixel_ind(1):replace_pixel_ind(1)+patch_size-1,...
        replace_pixel_ind(2):replace_pixel_ind(2)+patch_size-1,:);

    out_image_rgb(orig_pixel_index(1):orig_pixel_index(1)+patch_size-1,...
        orig_pixel_index(2):orig_pixel_index(2)+patch_size-1,:) = ...
        (cur_patch + replace_patch)/2;
end
out_image = lab2rgb(out_image_rgb);

end

