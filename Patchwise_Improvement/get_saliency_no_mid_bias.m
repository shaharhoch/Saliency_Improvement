function [ out_sal ] = get_saliency_no_mid_bias( in_img )
SALIENCY_ALGO_PATH = '../../PCA_Saliency_CVPR2013';

current_dir = cd; 
cd(SALIENCY_ALGO_PATH);
out_sal = PCA_Saliency_no_mid_bias(in_img);
cd(current_dir)

end

