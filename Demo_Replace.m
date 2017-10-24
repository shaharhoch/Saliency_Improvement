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
PATCH_SIZE = [7, 5, 3];

% Get image
image_path = fullfile(INPUT_FOLDER, INPUT_IMAGE);
in_image = im2double(imread(image_path));
in_image = rgb2lab(in_image);

% Get object mask
[~,im_name,im_ext] = fileparts(image_path);
mask_name = [im_name, '_Mask', im_ext];
mask_path = fullfile(MASK_FOLDER, mask_name);
object_mask = logical(imread(mask_path));

decimated_paths = cell(size(DECIMATION_RATIO));
saliency_paths = cell(size(DECIMATION_RATIO));
for patch_size_idx=1:length(PATCH_SIZE)
    patch_size = PATCH_SIZE(patch_size_idx);
    for decimination_idx=1:length(DECIMATION_RATIO)
        decimination_ratio = DECIMATION_RATIO(decimination_idx);
        out_image = in_image; 

        %Get replaced and left patches index
        [replaced_patches_idx, left_patches_idx] = getPatchesIdx(in_image,...
            object_mask, patch_size, decimination_ratio );

        for ind=1:length(replaced_patches_idx)
            closest_patch_ind = getClosesPatchInd(replaced_patches_idx(ind),...
                left_patches_idx, in_image, patch_size);

            %Replace the patch
            orig_pixel_index = patchIndToPixelInd(replaced_patches_idx(ind), patch_size, out_image);
            replace_pixel_ind = patchIndToPixelInd(closest_patch_ind, patch_size, out_image);

            out_image(orig_pixel_index(1):orig_pixel_index(1)+patch_size-1,...
                orig_pixel_index(2):orig_pixel_index(2)+patch_size-1,:) = ...
                in_image(replace_pixel_ind(1):replace_pixel_ind(1)+patch_size-1,...
                replace_pixel_ind(2):replace_pixel_ind(2)+patch_size-1,:);
        end
        out_image = lab2rgb(out_image);

        %Get image saliency
        current_dir = cd; 
        cd(SALIENCY_ALGO_PATH);
        saliency = PCA_Saliency(out_image);
        cd(current_dir)

        %Save the shifted image
        decimated_name = sprintf('%dx%d_%s_Decimated_%.1f', patch_size,...
            patch_size, im_name, decimination_ratio);
        decimated_path = fullfile(OUTPUT_FOLDER, [decimated_name, '.png']);
        decimated_paths{decimination_idx} = decimated_path;
        imwrite(out_image, decimated_path);

        %Save the saliency image
        saliency_name = sprintf('%dx%d_%s_Saliency_%.1f', patch_size,...
            patch_size, im_name, decimination_ratio);
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
        
        %Get saliency score
        saliency_score = getSaliencyScore(saliency, object_mask);

        subplot(2, length(DECIMATION_RATIO), ind+length(DECIMATION_RATIO));
        imshow(saliency)
        title(sprintf('Saliency Score=%.3f', saliency_score))
    end

    %Save the figure
    set(gcf,'Position',[0 0 1920 430]);
    
    fig_name = sprintf('%dx%d_%s_Decimation_Exam', patch_size, patch_size, im_name);
    savefig([OUTPUT_FOLDER, fig_name])
    export_fig(gcf, [OUTPUT_FOLDER, fig_name, '.png'])

    close all
end
