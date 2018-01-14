clc 
close all 
clear all

addpath('..\Export_Fig');

%Get info functions
INFO_FUNCS = {};
info_func.title = 'Saliency';
info_func.func = @(im_in) getSaliency(im_in);
INFO_FUNCS{1} = info_func;

info_func.title = 'Image Histograms';
info_func.func = @(im_in) getHistograms(im_in);
INFO_FUNCS{2} = info_func;

info_func.title = 'Image Gradient Map';
info_func.func = @(im_in) getGradients(im_in);
INFO_FUNCS{3} = info_func;

%Input/Output
INPUT_IMAGE = 'Car_Ad.jpg'; 
IMPROVED_INPUT_IMAGE = 'Car_Ad_Improved.jpg';
%INPUT_IMAGE = 'Perfume.jpg';

INPUT_FOLDER = '../Ad_Images/Images/';
MASK_FOLDER = '../Ad_Images/Masks/';
OUTPUT_FOLDER = '../Ad_Images/Out/';

% Get image
image_path = fullfile(INPUT_FOLDER, INPUT_IMAGE);
in_image = im2double(imread(image_path));

% Get image
image_path_improved = fullfile(INPUT_FOLDER, IMPROVED_INPUT_IMAGE);
in_image_improved = im2double(imread(image_path_improved));

% Get object mask
[~,im_name,im_ext] = fileparts(image_path);
mask_name = [im_name, '_Mask', im_ext];
mask_path = fullfile(MASK_FOLDER, mask_name);
object_mask = logical(imread(mask_path));

%Plot the info
figure(1);
set(gcf,'Visible','Off')
subplot(2, length(INFO_FUNCS)+1, 1)
imshow(in_image)
title('Original Image')

subplot(2, length(INFO_FUNCS)+1, length(INFO_FUNCS)+2)
imshow(in_image_improved)
title('Saliency Improved Image')

for ind=1:length(INFO_FUNCS)
    info_im_orig = INFO_FUNCS{ind}.func(in_image);
    info_im_improved = INFO_FUNCS{ind}.func(in_image_improved);
    
    figure(1)
    set(gcf,'Visible','Off')
    subplot(2, length(INFO_FUNCS)+1, ind+1)
    imshow(info_im_orig)
    title(INFO_FUNCS{ind}.title)
    
    subplot(2, length(INFO_FUNCS)+1, length(INFO_FUNCS)+2+ind)
    imshow(info_im_improved)
    title(INFO_FUNCS{ind}.title)
end

%Save the figure
set(gcf,'Position',[0 0 1920 400]);

fig_name = 'Saliency_Info_exam';
savefig([OUTPUT_FOLDER, fig_name])
export_fig(gcf, [OUTPUT_FOLDER, fig_name, '.jpg'],'-r400')
set(gcf,'Visible','On')

function [ out_sal ] = getHistograms( in_img )
    %Split into RGB Channels
    Red = in_img(:,:,1);
    Green = in_img(:,:,2);
    Blue = in_img(:,:,3);

    %Get histValues for each channel
    [yRed, x1] = imhist(Red);
    [yGreen, x2] = imhist(Green);
    [yBlue, x3] = imhist(Blue);
    
    % Normalize histograms
    yRed = yRed/sum(yRed);
    yGreen = yGreen/sum(yGreen);
    yBlue = yBlue/sum(yBlue);

    %Plot them together in one plot
    figure('visible','off')
    plot(x1, yRed, 'Red', x2, yGreen, 'Green', x3, yBlue, 'Blue');
    
    F = getframe(gcf);
    [out_sal, ~] = frame2im(F);
end

function [ out_grad ] = getGradients( in_img )
    im_gray = rgb2gray(in_img);
    [Gmag,Gdir] = imgradient(im_gray);
    
    Gx = Gmag.*cos(Gdir);
    Gy = Gmag.*sin(Gdir);
    
    %Plot them together in one plot
    figure('visible','off')
    quiver(Gx(end:-1:1,:), Gy(end:-1:1,:));
    axis tight
    set(gca, 'visible', 'off')
    
    F = getframe(gcf);
    [out_grad, ~] = frame2im(F);
end