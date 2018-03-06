clear all
close all
clc 
set(0,'DefaultFigureVisible','off');

addpath('..\Export_Fig');

%%% PARAMETERS DEFINITIONS %%%
OBJECT_MASK__TH = 230; 

IMG_REPLACE_OPT__ALPHA = 1500; 
IMG_REPLACE_OPT___REQ_IMPROVEMENT = 0.1;
IMG_REPLACE_OPT__PATCH_SIZE = 5; 

%%% PARAMETERS DEFINITIONS %%%
OUTPUT_FOLDER = '../Ad_Images/Out_Replace_Opt/';
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

avg_saliency_replace_opt = 0;
avg_saliency_orig = 0;
saliancy_score_diff = zeros(1, length(IN_IMAGES)); 
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
    
    % Get image after optimized replacement 
    [out_img_replace_opt, saliency_score_replace_opt] = image_patch_replace_opt( in_img, obj_mask,...
        IMG_REPLACE_OPT__ALPHA, IMG_REPLACE_OPT___REQ_IMPROVEMENT,...
        IMG_REPLACE_OPT__PATCH_SIZE);
    avg_saliency_replace_opt = avg_saliency_replace_opt + saliency_score_replace_opt;
    saliancy_score_diff(ind) = saliency_score_replace_opt - saliency_score_orig;
    
    out_replace_opt_name = sprintf('%s_replace_opt.png', in_img_name);
    out_replace_opt_path = fullfile(OUTPUT_FOLDER, out_replace_opt_name);
    imwrite(out_img_replace_opt, out_replace_opt_path);
end
avg_saliency_orig = avg_saliency_orig/length(IN_IMAGES);
avg_saliency_replace_opt = avg_saliency_replace_opt/length(IN_IMAGES);

% Get results bar graph 
names = {'Original', 'Optimized Patch Replacement'};
results = [avg_saliency_orig avg_saliency_replace_opt];
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

% Get saliency histogram
histogram(saliancy_score_diff)
title('Saliency Differnece Histogram'); 

%Save the figure
set(gcf,'Position',[0 0 700 900]);

fig_name = 'Saliency_Diff';
savefig([OUTPUT_FOLDER, fig_name])
export_fig(gcf, [OUTPUT_FOLDER, fig_name, '.png'])