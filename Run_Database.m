clear all
close all
clc 
set(0,'DefaultFigureVisible','off');

addpath('..\Export_Fig');

%%% PARAMETERS DEFINITIONS %%%
OBJECT_MASK__TH = 230; 

IMG_REPLACE__DECIMINATION_RATIO = 0.7; 
IMG_REPLACE__PATCH_SIZE = 7;

IMG_SHIFT__SHIFT_RATIO = 0.3;%0.2;

IMG_GAMMA__ALPHA = 0.1; 
IMG_GAMMA__GAMMA_0 = 1.5; 
IMG_GAMMA__GAMMA_OBJ = 0.6; 
IMG_GAMMA__FILTERED_SIZE = 10;
%%% PARAMETERS DEFINITIONS %%%
OUTPUT_FOLDER = '../Ad_Images/Out/';
IN_IMAGES_FOLDER = '../Ad_Images/IN_DB/';
IN_MASKS_FOLDER = '../Ad_Images/IN_DB/masks';

%Get DB list
IN_IMAGES = {};
IN_MASKS = {}; 

files_in_dir = dir(IN_IMAGES_FOLDER); 
for ind=1:length(files_in_dir)
    cur_file = files_in_dir(ind).name; 
    
    if(contains(cur_file,'jpg') == false)
        continue;
    end
    IN_IMAGES{end+1} = fullfile(IN_IMAGES_FOLDER, cur_file);
    
    [~, cur_file_name] = fileparts(cur_file);
    mask_file_name = sprintf('%s_mask.jpg', cur_file_name);
    IN_MASKS{end+1} = fullfile(IN_MASKS_FOLDER, mask_file_name);
end 
assert(isempty(IN_IMAGES) == false);
assert(length(IN_IMAGES) == length(IN_MASKS));

avg_saliency_orig = 0; 
avg_saliency_replace = 0; 
avg_saliency_shift = 0; 
avg_saliency_gamma = 0; 
avg_saliency_all = 0; 
saliancy_score_diff = []; 
for ind=1:length(IN_IMAGES)
    [~, in_img_name] = fileparts(IN_IMAGES{ind});
    in_img = im2double(imread(IN_IMAGES{ind}));
    
    if(ndims(in_img) ~= 3)
        continue; 
    end
    
    %Get mask
    obj_mask_tmp = im2uint8(imread(IN_MASKS{ind}));
    obj_mask = zeros(size(obj_mask_tmp));
    obj_mask(obj_mask_tmp > OBJECT_MASK__TH) = 1; 
    obj_mask = imresize(obj_mask,[size(in_img,1), size(in_img,2)]);
    obj_mask = logical(obj_mask);
    obj_mask = imfill(obj_mask, 'holes');
       
    %Save original image
    orig_name = sprintf('%s_orig.png', in_img_name);
    orig_path = fullfile(OUTPUT_FOLDER, orig_name);
    imwrite(in_img, orig_path);
    
    %Save mask
    mask_name = sprintf('%s_mask.png', in_img_name);
    mask_path = fullfile(OUTPUT_FOLDER, mask_name);
    imwrite(obj_mask, mask_path);
    
    % Get original saliency score
    saliency_orig = getSaliency(in_img);
    saliency_score_orig = getSaliencyScore(saliency_orig, obj_mask);
    avg_saliency_orig = avg_saliency_orig + saliency_score_orig;
    
    %Save original saliency
    imshow(saliency_orig);
    title(sprintf('Saliency Score = %f', saliency_score_orig));
    orig_saliency_name = sprintf('%s_orig_saliency.png', in_img_name);
    orig_saliency_path = fullfile(OUTPUT_FOLDER, orig_saliency_name);
    export_fig(gcf, orig_saliency_path);
    close;
    
    % Get image after replacement
    out_img_replace = image_replace_saliency(in_img, obj_mask,...
        IMG_REPLACE__DECIMINATION_RATIO, IMG_REPLACE__PATCH_SIZE );
    saliency_replace = getSaliency(out_img_replace);
    saliency_score_replace = getSaliencyScore(saliency_replace, obj_mask);
    avg_saliency_replace = avg_saliency_replace + saliency_score_replace;
    
    out_replace_name = sprintf('%s_replace.png', in_img_name);
    out_replace_path = fullfile(OUTPUT_FOLDER, out_replace_name);
    imwrite(out_img_replace, out_replace_path);
    
    % Get image after shift
    [out_img_shift, shifted_mask] = image_shift_saliency(in_img, obj_mask, IMG_SHIFT__SHIFT_RATIO);
    saliency_shift = getSaliency(out_img_shift);
    saliency_score_shift = getSaliencyScore(saliency_shift, shifted_mask);
    avg_saliency_shift = avg_saliency_shift + saliency_score_shift;
    
    out_shift_name = sprintf('%s_shift.png', in_img_name);
    out_shift_path = fullfile(OUTPUT_FOLDER, out_shift_name);
    imwrite(out_img_shift, out_shift_path);
    
    % Get image after gamma correction
    out_img_gamma = image_gamma_correction_saliency(in_img, obj_mask,...
        IMG_GAMMA__ALPHA, IMG_GAMMA__GAMMA_0, IMG_GAMMA__GAMMA_OBJ,...
        IMG_GAMMA__FILTERED_SIZE);
    saliency_gamma = getSaliency(out_img_gamma);
    saliency_score_gamma = getSaliencyScore(saliency_gamma, obj_mask);
    avg_saliency_gamma = avg_saliency_gamma + saliency_score_gamma;
    saliancy_score_diff = [saliancy_score_diff, (saliency_score_gamma-saliency_score_orig)];
    
    out_gamma_name = sprintf('%s_gamma.png', in_img_name);
    out_gamma_path = fullfile(OUTPUT_FOLDER, out_gamma_name);
    imwrite(out_img_gamma, out_gamma_path);
    
    imshow(saliency_gamma);
    title(sprintf('Saliency Score = %f', saliency_score_gamma));
    out_gamma_saliency_name = sprintf('%s_gamma_saliency.png', in_img_name);
    out_gamma_saliency_path = fullfile(OUTPUT_FOLDER, out_gamma_saliency_name);
    export_fig(gcf, out_gamma_saliency_path);
    close;
    
    % Get image after all
%     out_img_replace1 = out_img_replace;
%     out_img_gamma1 = image_gamma_correction_saliency(out_img_replace1, obj_mask,...
%         IMG_GAMMA__ALPHA, IMG_GAMMA__GAMMA_0, IMG_GAMMA__GAMMA_OBJ,...
%         IMG_GAMMA__FILTERED_SIZE);
%     [out_img_all, shifted_mask_all] = image_shift_saliency(out_img_gamma1,...
%         obj_mask, IMG_SHIFT__SHIFT_RATIO);
%     saliency_all = getSaliency(out_img_all);
%     saliency_score_all = getSaliencyScore(saliency_all, shifted_mask_all);
%     avg_saliency_all = avg_saliency_all + saliency_score_all;
%     
%     out_all_name = sprintf('%s_all.png', in_img_name);
%     out_all_path = fullfile(OUTPUT_FOLDER, out_all_name);
%     imwrite(out_img_all, out_all_path);
end

avg_saliency_orig = avg_saliency_orig/length(IN_IMAGES);
avg_saliency_replace = avg_saliency_replace/length(IN_IMAGES);
avg_saliency_shift = avg_saliency_shift/length(IN_IMAGES);
avg_saliency_gamma = avg_saliency_gamma/length(IN_IMAGES);
% avg_saliency_all = avg_saliency_all/length(IN_IMAGES);


% Get results bar graph 
% names = {'Original', 'Patch Replacement', 'Object Shift',...
%     'Gamma Correction', 'All'};
% results = [avg_saliency_orig avg_saliency_replace avg_saliency_shift...
%     avg_saliency_gamma avg_saliency_all];
names = {'Original', 'Patch Replacement', 'Object Shift',...
    'Gamma Correction'};
results = [avg_saliency_orig avg_saliency_replace avg_saliency_shift...
    avg_saliency_gamma];
bar(results)
set(gca,'xticklabel',names)
title('Average Saliency Score')
text(1:length(results),results,num2str(results',5),'vert','bottom','horiz','center'); 
box off

%Save the figure
set(gcf,'Position',[0 0 700 900]);

fig_name = 'Database_Run_Resutls';
savefig([OUTPUT_FOLDER, fig_name])
export_fig(gcf, [OUTPUT_FOLDER, fig_name, '.png'])

% Get gamma saliency histogram
histogram(saliancy_score_diff)
title('gamma_saliency-original_saliency'); 

%Save the figure
set(gcf,'Position',[0 0 700 900]);

fig_name = 'Saliency_Diff';
savefig([OUTPUT_FOLDER, fig_name])
export_fig(gcf, [OUTPUT_FOLDER, fig_name, '.png'])