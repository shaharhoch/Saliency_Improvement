function [ out_image ] = image_gamma_correction_saliency( in_img, object_mask, alpha, gamma_0, gamma_obj, filter_size )

%Constant mask for object and changing for background
%Alpah is a parameter that controls how quickly the gamma mask grows
%as the distance for the center of mass grows
gamma_mask = ones(size(object_mask));
gamma_mask(object_mask == 1) = gamma_obj; 

mask_center_of_mass = findMaskCenterOfMass(object_mask);
mask_area = sum(object_mask(:));
avg_radius = sqrt(mask_area/pi());
for i=1:size(gamma_mask, 1)
    for j=1:size(gamma_mask, 2)
        if(object_mask(i,j) == 1)
            continue; 
        end
        
        dist_center = sqrt((i-mask_center_of_mass(1))^2+(j-mask_center_of_mass(2))^2);
        gamma_mask(i,j) = gamma_0+alpha*((dist_center/avg_radius)-1); 
        
        if(gamma_mask(i,j) < 1)
            gamma_mask(i,j) = 1;
        end
    end
end
% Smooth the gamma mask 
avg_filter = ones(filter_size,filter_size) / (filter_size^2);
gamma_mask = imfilter(gamma_mask, avg_filter, 'circular');

%Correct the image
out_image = gamma_correction_mask(in_img, gamma_mask);
end

function [ center_of_mass ] = findMaskCenterOfMass(image_mask)
[row_mask, col_mask] = find(image_mask); 

center_of_mass = [mean(row_mask), mean(col_mask)];
center_of_mass = round(center_of_mass);
end

