clear all
close all
clc 

SALIENCY_ALGO_PATH = '../PCA_Saliency_CVPR2013';

%Add to path
addpath('StructCompletion\external\vlfeat-0.9.20\toolbox\mex\mexw64')
addpath('..\Export_Fig');

%Input/Output
INPUT_IMAGE = 'Car_Ad.jpg'; 
INPUT_IMAGE = 'Perfume.jpg';

INPUT_FOLDER = '../Ad_Images/Images/';
MASK_FOLDER = '../Ad_Images/Masks/';
OUTPUT_FOLDER = '../Ad_Images/Out/';

%The ratio of pixels that would be replaced
DECIMATION_RATIO = [0, 0.1, 0.2, 0.5, 0.7 0.9999999];
% Make sure 1 is not in the decimination ratio, because you can't replace
% all the pixels. 
assert(~ismember(1, DECIMATION_RATIO))

% Get image
image_path = fullfile(INPUT_FOLDER, INPUT_IMAGE);
in_image = im2double(imread(image_path));

% Get object mask
[~,im_name,im_ext] = fileparts(image_path);
mask_name = [im_name, '_Mask', im_ext];
mask_path = fullfile(MASK_FOLDER, mask_name);
object_mask = logical(imread(mask_path));

% Segment the image
%segments = vl_slic(single(rgb2lab(in_image)), 16, 300,'MinRegionSize',16);
segments = vl_slic(single(rgb2lab(in_image)), 40, 10, 'MinRegionSize',30);

% Make the masked part of the image be segment 0, so we won't touch it 
segments = segments + 1;
segments = segments .* uint32(abs(int32(object_mask)-1));

% Go over all the segments
decimated_paths = cell(size(DECIMATION_RATIO));
saliency_paths = cell(size(DECIMATION_RATIO));
for decimination_idx=1:length(DECIMATION_RATIO)
    decimination_ratio = DECIMATION_RATIO(decimination_idx);
    out_image = in_image; 
    for segment_idx=1:max(segments(:))
        segment_index = find(segments == segment_idx);
        number_decimated = floor(decimination_ratio*length(segment_index));
        
        replaced_idx = datasample(segment_index, number_decimated, 'Replace', false);
        left_idx = setdiff(segment_index, replaced_idx);
        
        % Find the replacement pixels index
        replacement_idx = datasample(left_idx, number_decimated, 'Replace', true);
        
        % Replace the pixels
        for ind=1:number_decimated
            [row_replaced, col_replaced] = ind2sub(size(segments), replaced_idx(ind));
            [row_replacement, col_replacement] = ind2sub(size(segments), replacement_idx(ind));
            out_image(row_replaced, col_replaced, :) = in_image(row_replacement, col_replacement, :);
        end
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
    im_title = sprintf('Deciminated Image Ratio: %.1f', decimination_ratio);
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
