clear all
close all
clc
fclose('all');

dbstop if error;
set(0,'DefaultFigureVisible','off');

addpath('..\..\Export_Fig');

run('options_def')

OBJECT_MASK__TH = 200;
OUTPUT_FOLDER = '../../Ad_Images/Out/';
IN_IMAGES_FOLDER = '../../Ad_Images/IN_DB/';
IN_MASKS_FOLDER = '../../Ad_Images/IN_DB/masks';

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

% Init log file 
log_path = fullfile(OUTPUT_FOLDER, 'Saliency.log');
log_fid = fopen(log_path, 'w');

% [orig, ours, roey]
saliency_arr = []; 
saliency_score_roey = 0; %Temporary.
for ind=1:length(IN_IMAGES)    
    [~, in_img_name] = fileparts(IN_IMAGES{ind});
    in_img = im2double(imread(IN_IMAGES{ind}));
    
    % Make sure this is an RGB image
    if(ndims(in_img) ~= 3)
        continue; 
    end
    
    %Get mask
    obj_mask_tmp = im2uint8(imread(IN_MASKS{ind}));
    obj_mask_tmp = imresize(obj_mask_tmp,[size(in_img,1), size(in_img,2)]);
    obj_mask = zeros(size(obj_mask_tmp));
    obj_mask(obj_mask_tmp > OBJECT_MASK__TH) = 1; 
    obj_mask = logical(obj_mask);
        
    %Save mask
    mask_name = sprintf('%s_mask.png', in_img_name);
    mask_path = fullfile(OUTPUT_FOLDER, mask_name);
    imwrite(obj_mask, mask_path);
    
    % Run algorithm 
    [ out_ours, saliency_ours, saliency_orig, num_iter ] = Patchwise_Improve( in_img, obj_mask );
    saliency_score_orig = getSaliencyScore(saliency_orig, obj_mask);
    saliency_score_ours = getSaliencyScore(saliency_ours, obj_mask);
    
    % Run Roey's algorithm
%     out_roey = improve_saliency_roey(IN_IMAGES{ind}, IN_MASKS{ind});
%     saliency_roey = getSaliency(out_roey);
%     saliency_score_roey = getSaliencyScore(saliency_roey, obj_mask);
    
    %Save original image and saliency
    orig_name = sprintf('%s_img_orig.png', in_img_name);
    orig_path = fullfile(OUTPUT_FOLDER, orig_name);
    imwrite(in_img, orig_path);
    
    imshow(saliency_orig);
    title(sprintf('Saliency Score = %f', saliency_score_orig));
    our_saliency_name = sprintf('%s_saliency_orig.png', in_img_name);
    our_saliency_path = fullfile(OUTPUT_FOLDER, our_saliency_name);
    export_fig(gcf, our_saliency_path);
    close;
    
    % Save our image and saliency    
    out_ours_name = sprintf('%s_img_ours.png', in_img_name);
    out_ours_path = fullfile(OUTPUT_FOLDER, out_ours_name);
    imwrite(out_ours, out_ours_path);
    
    imshow(saliency_ours);
    title(sprintf('Saliency Score = %f', saliency_score_ours));
    our_saliency_name = sprintf('%s_saliency_our.png', in_img_name);
    our_saliency_path = fullfile(OUTPUT_FOLDER, our_saliency_name);
    export_fig(gcf, our_saliency_path);
    close;
    
%     %Save Roey image and saliency
%     roey_name = sprintf('%s_img_roey.png', in_img_name);
%     roey_path = fullfile(OUTPUT_FOLDER, roey_name);
%     imwrite(in_img, roey_path);
    
%     imshow(saliency_roey);
%     title(sprintf('Saliency Score = %f', saliency_score_roey));
%     roey_saliency_name = sprintf('%s_saliency_roey.png', in_img_name);
%     roey_saliency_path = fullfile(OUTPUT_FOLDER, roey_saliency_name);
%     export_fig(gcf, roey_saliency_path);
%     close;
%     
    % Update saliency array
    saliency_arr = [saliency_arr; [saliency_score_orig, saliency_score_ours, saliency_score_roey]];
    
    % Update saliency in output file
    fprintf(log_fid, '%s: Original = %3f, Ours = %3f, Roey = %3f, Num of iterations: %d\n',...
        in_img_name, saliency_score_orig, saliency_score_ours, saliency_score_roey, num_iter);
end

fclose(log_fid); 

% Save the saliency array
log_saliency = fullfile(OUTPUT_FOLDER, 'Saliency_Scores');
save(log_saliency, 'saliency_arr');

