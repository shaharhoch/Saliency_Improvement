function [ out_img, cur_saliancy_score, num_iter ] = image_patch_replace_opt( in_img, object_mask, alpha, req_improvement, patch_size )
% This function performs an image patch replacement optimization to 
% improve the saliency of an image.
% in_img- the input image.
% alpha- the ratio of importance between keeping the image as is and the
%        saliency improvement.
% req_improvement- the required improvement in saliency score.
% K- number of nearset neigbours.

% Note about patches: 
% The way we treat patches is that each patch has an index, starting from
% the top left corener, moving right and then down. This way each patch can
% be refered to using only one index.

%Remove from the image the parts that are not fully divided into patches
new_img_size = floor(size(in_img)/patch_size)*patch_size;
in_img = in_img(1:new_img_size(1), 1:new_img_size(2), :);
object_mask = object_mask(1:new_img_size(1), 1:new_img_size(2));

%Make sure the image is double
in_img = im2double(in_img); 

saliency = getSaliency(in_img);
orig_saliency_score = getSaliencyScore(saliency, object_mask);
saliency_improvement = 0; 
cur_saliancy_score = orig_saliency_score;
num_iter = 0; 
out_img = in_img;
while(saliency_improvement < req_improvement)
    out_img_tmp = single_iteration(in_img, object_mask, alpha, patch_size, saliency);
    num_iter = num_iter + 1;
    
    saliency = getSaliency(out_img_tmp);
    new_saliancy_score = getSaliencyScore(saliency, object_mask);
    if(new_saliancy_score <= cur_saliancy_score)
        break; 
    end
    
    out_img = out_img_tmp;
    cur_saliancy_score = new_saliancy_score;
    saliency_improvement = cur_saliancy_score - orig_saliency_score; 
end
    
end

function [ out_img ] = single_iteration( in_img, object_mask, alpha, patch_size, saliency )
%Get lab coordinates for the image
in_img_lab = rgb2lab(in_img);

%Get a matrix where each row is a patch and the columns are the lab values
%of the image.
img_size = size(in_img);
num_of_patches = (img_size(1)*img_size(2))/(patch_size^2);
lab_data = zeros([num_of_patches, 3*(patch_size^2)]);
for i=1:num_of_patches
    mtx_ind = patch_ind_to_mtx_ind(i, patch_size, in_img);
    lab_data(i,:) = reshape(in_img_lab(mtx_ind(1):(mtx_ind(1)+patch_size-1),...
        mtx_ind(2):(mtx_ind(2)+patch_size-1), :), 1, size(lab_data, 2));
end

% Get the extended data with the normalized saliency as the last column
lab_data_ext = [lab_data, zeros(size(lab_data,1), 1)];
lab_data_with_sal = lab_data_ext; 
for i=1:num_of_patches
    lab_data_with_sal(i, end) = alpha*get_patch_saliency(saliency, i, patch_size);
end
   
% Find KNN
nn_idx = knnsearch(lab_data_with_sal, lab_data_ext); 

% Replace patches
out_img = in_img;
for i=1:num_of_patches
    if(nn_idx(i) == i)
        continue
    end
    
    cur_patch_idx = patch_ind_to_mtx_ind(i, patch_size, in_img);
    if(object_mask(cur_patch_idx(1), cur_patch_idx(2)) == 1)
        %Don't replace pixels inside the object mask
        continue
    end
    
    rep_patch_idx = patch_ind_to_mtx_ind(nn_idx(i), patch_size, in_img);
    out_img(cur_patch_idx(1):(cur_patch_idx(1)+patch_size-1),...
        cur_patch_idx(2):(cur_patch_idx(2)+patch_size-1), :) = ...
        in_img(rep_patch_idx(1):(rep_patch_idx(1)+patch_size-1),...
        rep_patch_idx(2):(rep_patch_idx(2)+patch_size-1), :);
end

end


function [ mtx_ind ] = patch_ind_to_mtx_ind( patch_ind, patch_size, in_img )
    mtx_ind = zeros(1,2);
    mtx_ind(1) = 1+patch_size*floor(((patch_ind-1)*patch_size)/size(in_img,2));
    mtx_ind(2) = mod(((patch_ind-1)*patch_size), size(in_img,2)) + 1;
end

function [ patch_saliency ] = get_patch_saliency( saliency, patch_ind, patch_size )
    mtx_ind = patch_ind_to_mtx_ind(patch_ind, patch_size, saliency);
    
    patch = saliency(mtx_ind(1):(mtx_ind(1)+patch_size-1), mtx_ind(2):(mtx_ind(2)+patch_size-1));
    patch_saliency = sum(patch(:))/(patch_size^2);
end
