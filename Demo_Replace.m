clear all
close all
clc 

SALIENCY_ALGO_PATH = '../PCA_Saliency_CVPR2013';

%Add to path
addpath('StructCompletion\external\vlfeat-0.9.20\toolbox\mex\mexw64')
addpath('..\Export_Fig');

%Input/Output
INPUT_IMAGE = 'Car_Ad.jpg'; 
%INPUT_IMAGE = 'Perfume.jpg';

INPUT_FOLDER = '../Ad_Images/Images/';
MASK_FOLDER = '../Ad_Images/Masks/';
OUTPUT_FOLDER = '../Ad_Images/Out/';

%The ratio of pixels that would be replaced
DECIMATION_RATIO = [0, 0.1, 0.2, 0.5, 0.7 0.99];
% Make sure 1 is not in the decimination ratio, because you can't replace
% all the pixels. 
assert(~ismember(1, DECIMATION_RATIO))
PATCH_SIZE = 5;

% Get image
image_path = fullfile(INPUT_FOLDER, INPUT_IMAGE);
in_image = im2double(imread(image_path));

% Get object mask
[~,im_name,im_ext] = fileparts(image_path);
mask_name = [im_name, '_Mask', im_ext];
mask_path = fullfile(MASK_FOLDER, mask_name);
object_mask = logical(imread(mask_path));

decimated_paths = cell(size(DECIMATION_RATIO));
saliency_paths = cell(size(DECIMATION_RATIO));
for decimination_idx=1:length(DECIMATION_RATIO)
    decimination_ratio = DECIMATION_RATIO(decimination_idx);
    out_image = in_image; 
    
    %Get replaced and left patches index
    [replaced_patches_idx, left_patches_idx] = getPatchesIdx(in_image,...
        object_mask, PATCH_SIZE, decimination_ratio );
    
    for ind=1:length(replaced_patches_idx)
        closest_patch_ind = getClosesPatchInd(replaced_patches_idx(ind),...
            left_patches_idx, in_image, PATCH_SIZE);
        
        %Replace the patch
        orig_pixel_index = patchIndToPixelInd(replaced_patches_idx(ind), PATCH_SIZE, out_image);
        replace_pixel_ind = patchIndToPixelInd(closest_patch_ind, PATCH_SIZE, out_image);
        
        out_image(orig_pixel_index(1):orig_pixel_index(1)+PATCH_SIZE-1,...
            orig_pixel_index(2):orig_pixel_index(2)+PATCH_SIZE-1,:) = ...
            in_image(replace_pixel_ind(1):replace_pixel_ind(1)+PATCH_SIZE-1,...
            replace_pixel_ind(2):replace_pixel_ind(2)+PATCH_SIZE-1,:);
    end
    
    %Get image saliency
    current_dir = cd; 
    cd(SALIENCY_ALGO_PATH);
    saliency = PCA_Saliency(out_image);
    cd(current_dir)
    
    %Save the shifted image
    decimated_name = sprintf('%s_Decimated_%.1f', im_name, decimination_ratio, im_ext);
    decimated_path = fullfile(OUTPUT_FOLDER, [decimated_name, '.png']);
    decimated_paths{decimination_idx} = decimated_path;
    imwrite(out_image, decimated_path);
    
    %Save the saliency image
    saliency_name = sprintf('%s_Saliency_%.1f', im_name, decimination_ratio, im_ext);
    saliency_path = fullfile(OUTPUT_FOLDER, [saliency_name, '.png']);
    imwrite(saliency, saliency_path);
    saliency_paths{decimination_idx} = saliency_path;
end

% Plot the figure
figure;
for ind=1:length(decimated_paths)
    decimination_ratio = DECIMATION_RATIO(ind);
    image = imread(decimated_paths{ind});
    saliency = imread(saliency_paths{ind});
    
    hold on;
    subplot(2, length(DECIMATION_RATIO), ind);
    imshow(image)
    im_title = sprintf('Replacement Ratio: %.1f', decimination_ratio);
    title(im_title)
    
    subplot(2, length(DECIMATION_RATIO), ind+length(DECIMATION_RATIO));
    imshow(saliency)
    title('Saliency')
end

%Save the figure
set(gcf,'Position',[0 0 1920 430]);

savefig([OUTPUT_FOLDER, im_name, '_Decimation_Exam'])
export_fig(gcf, [OUTPUT_FOLDER, im_name, '_Decimation_Exam.png'])

close all
