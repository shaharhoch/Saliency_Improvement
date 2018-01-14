clear all
close all
clc 

SALIENCY_ALGO_PATH = '../PCA_Saliency_CVPR2013';
addpath('criminisi_inpainting');
addpath(SALIENCY_ALGO_PATH);
addpath('..\Export_Fig');
addpath(cd)

%%%%% DEFINES %%%%%
%Extrapolation
%Type can be either 'STRUCT_COMPLETION' or 'CRIMMINISI'
EXTRAPOLATION_TYPE = 'STRUCT_COMPLETION';
%Crimminisi
PATCH_SIZE = 9;

%Interpolation
INTERPOLATION_RANDOMNESS = 0;%9;
INTERPOLATION_RATIO = 0;

%Input/Output
%INPUT_IMAGE = 'Car_Ad.jpg'; 
INPUT_IMAGE = 'Perfume.jpg';

INPUT_FOLDER = '../Ad_Images/Images/';
MASK_FOLDER = '../Ad_Images/Masks/';
OUTPUT_FOLDER = '../Ad_Images/Out/';
SHIFT_RATIOS = [0,0.2,0.5,0.8,1]; 

% Get image
image_path = fullfile(INPUT_FOLDER, INPUT_IMAGE);
in_image = im2double(imread(image_path));

% Get object mask
[~,im_name,ext] = fileparts(image_path);
mask_name = [im_name, '_Mask', ext];
mask_path = fullfile(MASK_FOLDER, mask_name);
object_mask = logical(imread(mask_path));

shift_paths = cell(length(SHIFT_RATIOS));
saliency_paths = cell(length(SHIFT_RATIOS));
shifted_mask_paths = cell(length(SHIFT_RATIOS));
for ind=1:length(SHIFT_RATIOS)
    shift_ratio = SHIFT_RATIOS(ind);
    
    %Shift the image
    [shifted_im, new_image_bounds, shifted_mask] = imageShift(in_image, object_mask, shift_ratio);
    
    %%%%%%%%%%% Fill the image using interpolation %%%%%%%%%%%
    %Get wanted interpolation boarders
    interpolation_boarders = new_image_bounds;
    interpolation_boarders(:,1) = new_image_bounds(:,1) -...
        round(INTERPOLATION_RATIO*[new_image_bounds(1,1)-1; new_image_bounds(2,1)-1]);
    interpolation_boarders(:,2) = new_image_bounds(:,2) +...
        round(INTERPOLATION_RATIO*[size(shifted_im,1)-new_image_bounds(1,2); size(shifted_im,2)-new_image_bounds(2,2)]);
            
    interpolated_image = interpolateImage(shifted_im, new_image_bounds, INTERPOLATION_RANDOMNESS, interpolation_boarders);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%% Fill the image using impainting %%%%%%%%%%%
    if(strcmp(EXTRAPOLATION_TYPE, 'STRUCT_COMPLETION'))
        % Move working directory to the struct completion directory
        save_cd = cd;
        cd('StructCompletion')
        
        %save tmp_image
        mask = ImageBoardersToMask(interpolated_image, interpolation_boarders);  
        tmp_image_name = 'tmp_img.png';
        imwrite(interpolated_image, ['data\',tmp_image_name], 'Alpha', mask)
        
        %Run the extrapolation
        filled_image = sc_complete(tmp_image_name);
        
        %Return working directory
        cd(save_cd)
    elseif(strcmp(EXTRAPOLATION_TYPE, 'CRIMMINISI'))
        filled_image = crimminisiExtrapolate(interpolated_image, interpolation_boarders, PATCH_SIZE);
    else
        error('Invalid EXTRAPOLATION_TYPE')
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    %Get image saliency
    current_dir = cd; 
    cd(SALIENCY_ALGO_PATH);
    saliency = PCA_Saliency(filled_image);
    cd(current_dir)
    
    %Save the shifted image
    shifted_name = sprintf('%s_Shifted_%.1f', im_name, shift_ratio, ext);
    shifted_path = fullfile(OUTPUT_FOLDER, [shifted_name, '.jpg']);
    shift_paths{ind} = shifted_path;
    imwrite(filled_image, shifted_path);
    
    %Save the saliency image
    saliency_name = sprintf('%s_Saliency_%.1f', im_name, shift_ratio, ext);
    saliency_path = fullfile(OUTPUT_FOLDER, [saliency_name, '.jpg']);
    imwrite(saliency, saliency_path);
    saliency_paths{ind} = saliency_path;
    
    %Save the shifted mask
    shifted_mask_name = sprintf('%s_Shifted_Mask_%.1f', im_name, shift_ratio, ext);
    shifted_mask_path = fullfile(OUTPUT_FOLDER, [shifted_mask_name, '.jpg']);
    imwrite(shifted_mask, shifted_mask_path);
    shifted_mask_paths{ind} = shifted_mask_path;
end

% Plot the figure
figure;
for ind=1:length(shift_paths)
    shift_ratio = SHIFT_RATIOS(ind);
    image = imread(shift_paths{ind});
    saliency = imread(saliency_paths{ind});
    shifted_mask = im2single(imread(shifted_mask_paths{ind}));
    
    saliency_score = getSaliencyScore(saliency, shifted_mask);
    
    hold on;
    subplot(2, length(SHIFT_RATIOS), ind);
    imshow(image)
    im_title = sprintf('Shifted Image Ratio: %.1f', shift_ratio);
    title(im_title)
    
    subplot(2, length(SHIFT_RATIOS), ind+length(SHIFT_RATIOS));
    imshow(saliency)
    title('Saliency')
    
    saliency_score_label = sprintf('Saliency Score: %.2f', saliency_score);
    xlabel(saliency_score_label)
end

%Save the figure
set(gcf,'Position',[0 0 1920 430]);

savefig([OUTPUT_FOLDER, im_name, '_Shifted_Exam'])
export_fig(gcf, [OUTPUT_FOLDER, im_name, '_Shifted_Exam.png'])

close all