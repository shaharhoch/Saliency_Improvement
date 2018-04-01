clear all
close all
clc

run('options_def')

in_img = imread('C:\Users\shaha\Dropbox\School\Project_Ayelet\Ad_Images\IN_DB\108.jpg');
in_mask = imread('C:\Users\shaha\Dropbox\School\Project_Ayelet\Ad_Images\IN_DB\masks\108_mask.jpg');
im_out = Patchwise_Improve(in_img, in_mask);

imshow(im_out)