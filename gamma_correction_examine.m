clear all
close all
clc 

addpath('..\Export_Fig');

SALIENCY_ALGO_PATH = '../PCA_Saliency_CVPR2013';

%Input/Output
%INPUT_IMAGE = 'Car_Ad.jpg'; 
INPUT_IMAGE = 'Perfume.jpg';

INPUT_FOLDER = '../Ad_Images/Images/';
MASK_FOLDER = '../Ad_Images/Masks/';
OUTPUT_FOLDER = '../Ad_Images/Out_Local/';

% Get image
image_path = fullfile(INPUT_FOLDER, INPUT_IMAGE);
in_image = im2double(imread(image_path));

% Get object mask
[~,im_name,im_ext] = fileparts(image_path);
mask_name = [im_name, '_Mask', im_ext];
mask_path = fullfile(MASK_FOLDER, mask_name);
object_mask = logical(imread(mask_path));

%Create Gamma masks
GAMMA_MASKS = {}; 

% No Change mask
gamma_mask = ones(size(object_mask));
GAMMA_MASKS{1} = gamma_mask;

% Constant mask for object and constant for background
gamma_mask = ones(size(object_mask));
gamma_mask(object_mask == 1) = 0.4; 
gamma_mask(object_mask == 0) = 2;
% Smooth the gamma mask 
avg_filter = ones(80,80) / (80^2);
gamma_mask = imfilter(gamma_mask, avg_filter, 'circular');
GAMMA_MASKS{2} = gamma_mask;

%Constant mask for object and changing for background
 %Alpah is a parameter that controls how quickly the gamma mask grows
 %as the distance for the center of mass grows
alpha = 0.7;
gamma_mask = ones(size(object_mask));
gamma_mask(object_mask == 1) = 0.4; 

mask_center_of_mass = findMaskCenterOfMass(object_mask);
mask_area = sum(object_mask(:));
avg_radius = sqrt(mask_area/pi());
for i=1:size(gamma_mask, 1)
    for j=1:size(gamma_mask, 2)
        if(object_mask(i,j) == 1)
            continue; 
        end
        
        dist_center = sqrt((i-mask_center_of_mass(1))^2+(j-mask_center_of_mass(2))^2);
        gamma_mask(i,j) = alpha*(dist_center/avg_radius); 
        
        if(gamma_mask(i,j) < 1)
            gamma_mask(i,j) = 1;
        end
    end
end
% Smooth the gamma mask 
avg_filter = ones(80,80) / (80^2);
gamma_mask = imfilter(gamma_mask, avg_filter, 'circular');
GAMMA_MASKS{3} = gamma_mask;

figure(1); 
for ind=1:length(GAMMA_MASKS)
    gamma_corrected_img = gamma_correction_mask(in_image, GAMMA_MASKS{ind}); 
    
    %Get image saliency
    current_dir = cd; 
    cd(SALIENCY_ALGO_PATH);
    saliency = PCA_Saliency(gamma_corrected_img);
    cd(current_dir)
    saliency = im2uint8(saliency);
    
    %Plot corrected image 
    subplot(length(GAMMA_MASKS), 3, (ind-1)*3+1); 
    imshow(gamma_corrected_img); 
    title('Gamma Corrected Image')
    
    %Plot saliency
    subplot(length(GAMMA_MASKS), 3, (ind-1)*3+2); 
    imshow(saliency);
    saliency_score = getSaliencyScore(saliency, object_mask);
    xlabel(sprintf('Saliency Score = %.2f', saliency_score));
    title('Saliency Map')
    
    %Plot Gamma Correction map
    subplot(length(GAMMA_MASKS), 3, (ind-1)*3+3); 
    imagesc(GAMMA_MASKS{ind})
    colorbar
    title('Gamma Correction Map')
end

%Save the figure
set(gcf,'Position',[0 0 1000 1080]);

timestamp = datestr(datetime('now'), 'dd_mm_yy__HH_MM_SS');
fig_name = sprintf('Gamma_Correction_Exam_%s', timestamp);
savefig([OUTPUT_FOLDER, fig_name])
export_fig(gcf, [OUTPUT_FOLDER, fig_name, '.png'])

function [ center_of_mass ] = findMaskCenterOfMass(image_mask)
[row_mask, col_mask] = find(image_mask); 

center_of_mass = [mean(row_mask), mean(col_mask)];
center_of_mass = round(center_of_mass);
end
    