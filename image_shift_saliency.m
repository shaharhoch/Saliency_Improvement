function [ out_image, shifted_mask ] = image_shift_saliency( in_img, object_mask, shift_ratio )

%Shift the image
[shifted_im, new_image_bounds, shifted_mask] = imageShift(in_img, object_mask, shift_ratio);
mask_img = ImageBoardersToMask(in_img, new_image_bounds);  

%%%%%%%%%%% Fill the image using image melding %%%%%%%%%%%
% Move working directory to the image melding directory
% save_cd = cd;
% cd('Image_Melding')
% 
% tic
% mask_hole = 1-mask_img;
% out_image = HoleFillingFunc( shifted_im, mask_hole );
% toc;
% 
% %Return to working directory
% cd(save_cd)


%%%%%%%%%%% Fill the image using impainting %%%%%%%%%%%
% Move working directory to the struct completion directory
save_cd = cd;
cd('C:\StructCompletion')

%Save temp image for struct completion
tmp_image_name = 'tmp_img.png';
imwrite(shifted_im, ['data\',tmp_image_name], 'Alpha', mask_img)

%Run the extrapolation
out_image = sc_complete(tmp_image_name);
close

%Return to working directory
cd(save_cd)
end

