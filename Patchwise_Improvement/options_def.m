global options

options.percentage_of_processed_pixels_in_interation_object = 100;
options.percentage_of_processed_pixels_in_interation_background = 50;
options.patch_size = 7;
options.max_num_iter = 5;
options.req_sal_score_improvement = 0.15;

%Alpha os the ratio in which you multiply a patch in pca space. If it is
%larger than 1 it increases saliency and if it is smaller than 1 it loweres
%saliency. Alpha needs to be a 1X3 vector, each value is for a different
%channel in L, a, b
options.mask_alpha = 1.15*[1, 1, 1]; 
options.background_alpha = 0.85*[1, 1, 1];

options.weighted_alpha = false; % Needs to be true/false
options.pca_dim_th_ration = 0.005;

% Weights of gradients for SPE, it needs to be 3 channels for  L, a, b
options.spe_grad_weight = [5, 2, 2];