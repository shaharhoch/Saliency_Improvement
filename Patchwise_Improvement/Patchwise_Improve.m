function [ out_im ] = Patchwise_Improve( in_img, in_mask )
global options

in_img = im2double(in_img);

% Fix mask 
MASK_TH = 0.8;
in_mask_double = im2double(in_mask);
in_mask_tmp = false(size(in_mask));
in_mask_tmp(in_mask_double > MASK_TH) = true; 
in_mask = in_mask_tmp;

in_mask = imresize(in_mask,[size(in_img, 1), size(in_img, 2)]);

% Get gradients
in_lab = rgb2lab(in_img); 
in_grad = cell(3,2);
for i=1:3
    [dx, dy] = sGradMex(single(in_lab(:,:,i))); 
    in_grad{i, 1} = dx;
    in_grad{i, 2} = dy; 
end

out_im = in_img;
for i=1:options.num_iter
    out_im = Patchwise_Improve_Iteration(out_im, in_mask, in_lab, in_grad);
end

end

function [ out_im ] = Patchwise_Improve_Iteration( in_img, in_mask, in_lab, in_grad )
% Perform a patch-wise improvement of the saliency of the object marked by
% the mask.
global options

[ max_sal_ind,  min_sal_ind ] = get_min_max_sal_ind(in_img, in_mask);

% Manipulate patcehs with the highest saliency in background
out_lab = in_lab;
for i=1:3
    out_lab(:, :, i) = modify_channel_patches(out_lab(:, :, i), max_sal_ind,...
        options.background_alpha(i));
end

% Manipulate patcehs with the lowest saliency in the object
for i=1:3
    out_lab(:, :, i) = modify_channel_patches(out_lab(:, :, i), min_sal_ind,...
        options.mask_alpha(i));
end

% Use Screened Poisson Equation (SPE) to make the image gradients look more
% natural
for i=1:3
    out_lab(:, :, i) = doPoissonCombination(out_lab(:, :, i),...
        in_grad{i,1}, in_grad{i,2}, options.spe_grad_weight(i));
end

out_im = lab2rgb(out_lab);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_min_max_sal_ind- Get the index of patches to modify.
% max_sal_ind- Indexes of most salient pixels in the background. 
% min_sal_ind- Indexes of least salient pixels in the object.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ max_sal_ind,  min_sal_ind ] = get_min_max_sal_ind( in_img, in_mask )
global options

% Get saliency 
saliency = get_saliency_no_mid_bias(in_img); 

num_of_processed_patches = round((options.percentage_of_processed_pixels_in_interation*numel(saliency))/100);
% Get indexes of lowest saliency pixels in the object
num_of_processed_patches_in_object = round((options.percentage_of_processed_pixels_in_interation*sum(in_mask(:)))/100); 
saliency_masked = double(saliency); 
saliency_masked(in_mask ~= 1) = 300;
[~, min_sal_ind] = sort(saliency_masked(:), 'ascend');
min_sal_ind = min_sal_ind(1:num_of_processed_patches_in_object);

% Get indexes of highest saliency pixels in the backgorund
num_of_processed_patches_in_background = num_of_processed_patches - num_of_processed_patches_in_object;
saliency_anti_masked = double(saliency); 
saliency_anti_masked(in_mask == 1) = -1;
[~, max_sal_ind] = sort(saliency_anti_masked(:), 'descend');
max_sal_ind = max_sal_ind(1:num_of_processed_patches_in_background);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% modify_channel_patches- Modify patches of a single channel.
% im_channel- One channel of an image.
% patches_idx- The patches to modify. 
% alpha- The parameter of change of the patches. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ channel_out ] = modify_channel_patches(im_channel, patches_idx, alpha)
% Get PCA of patches
global options 

patch_size = [options.patch_size, options.patch_size];
patches_mtx = im2patches(im_channel, patch_size);

mean_patch = mean(patches_mtx,1);
mean_mtx = repmat(mean_patch, [size(patches_mtx, 1), 1]);
[coeff,score] = pca(patches_mtx - mean_mtx);

for i=1:length(patches_idx)
    ind = patches_idx(i); 
    
    score(ind, :) = score(ind, :) * alpha;
end

% Go from patches to a normal image
patches_mtx_after_pca = (score * coeff') + mean_mtx;
channel_out = patches2im(patches_mtx_after_pca, patch_size, size(im_channel));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% patches2im- Converts patches representation to an image by 
% averaging all contributions from patches that touch it. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ out_img ] = patches2im(in_patches, patch_size, out_im_size)
assert(mod(patch_size(1), 2) == 1);
assert(mod(patch_size(2), 2) == 1);
assert(out_im_size(1)*out_im_size(2) == size(in_patches, 1));

out_img = zeros(out_im_size);
num_of_patches = zeros(out_im_size);

row_half_patch = (patch_size(1)-1)/2;
col_half_patch = (patch_size(2)-1)/2;

for ind=1:size(in_patches,1)
    [row, col] = ind_to_col_row(ind, out_im_size);
    
    start_row = max(row-row_half_patch, 1);
    end_row = min(row+row_half_patch, out_im_size(1));
    
    start_col = max(col-col_half_patch, 1);
    end_col = min(col+col_half_patch, out_im_size(2));
    
    patch = reshape(in_patches(ind, :), patch_size);
    
    patch_start_row = (row_half_patch + 1) - (row - start_row);
    patch_end_row = (row_half_patch + 1) + (end_row - row);
    
    patch_start_col = (col_half_patch + 1) - (col - start_col);
    patch_end_col = (col_half_patch + 1) + (end_col - col);
    
    out_img(start_row:end_row, start_col:end_col) = ...
        out_img(start_row:end_row, start_col:end_col) + ...
        patch(patch_start_row:patch_end_row, patch_start_col:patch_end_col);
    
    num_of_patches(start_row:end_row, start_col:end_col) = ...
        num_of_patches(start_row:end_row, start_col:end_col) + 1;
end
num_of_patches(num_of_patches == 0) = 1;

out_img = out_img./num_of_patches;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% im2patches- Converts an image (only one channel) to a patch
% representation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ out_patches ] = im2patches(in_img, patch_size)
assert(mod(patch_size(1), 2) == 1);
assert(mod(patch_size(2), 2) == 1);

row_half_patch = (patch_size(1)-1)/2;
col_half_patch = (patch_size(2)-1)/2;

out_patches = zeros(numel(in_img), patch_size(1)*patch_size(2));
in_img_size = size(in_img);

for row=1:size(in_img, 1)
    for col=1:size(in_img, 2)
        ind = row + (size(in_img, 1)*(col - 1));
        
        start_row = max(row-row_half_patch, 1);
        end_row = min(row+row_half_patch, in_img_size(1));

        start_col = max(col-col_half_patch, 1);
        end_col = min(col+col_half_patch, in_img_size(2));

        patch_start_row = (row_half_patch + 1) - (row - start_row);
        patch_end_row = (row_half_patch + 1) + (end_row - row);

        patch_start_col = (col_half_patch + 1) - (col - start_col);
        patch_end_col = (col_half_patch + 1) + (end_col - col);
        
        patch = zeros(patch_size) + in_img(row, col);
        patch(patch_start_row:patch_end_row, patch_start_col:patch_end_col) = ...
            in_img(start_row:end_row, start_col:end_col);
        
        out_patches(ind, :) = reshape(patch, [1, patch_size(1)*patch_size(2)]);
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ind_to_col_row- Converts a single index to 2 index of row 
% and column.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [row, col] = ind_to_col_row(ind, im_size)
row = mod(ind-1, im_size(1))+1;
col = ceil(ind/im_size(1));
end

