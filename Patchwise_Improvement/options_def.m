global options

options.percentage_of_processed_pixels_in_interation = 100;
options.patch_size = 7;

%Alpha os the ratio in which you multiply a patch in pca space. If it is
%larger than 1 it increases saliency and if it is smaller than 1 it loweres
%saliency. Alpha needs to be a 1X3 vector, each value is for a different
%channel in L, a, b
options.mask_alpha = 4*[1, 1, 1]; 
options.background_alpha = 0.2*[1, 1, 1];

options.num_iter = 4;

% Weights of gradients for SPE, it needs to be 3 channels for  L, a, b
options.spe_grad_weight = [3, 1, 1];