function [ image_out, new_image_bounds, shifted_mask ] = imageShift( image, image_mask, shift_ratio )
%This function shifts an image, so the desired object will be closer to the
%center.
% image- the original image
% image_mask- a mask of the image you want to center.
% shift_amount- a number between 0 and 1 that says how much to shift the
% object toward the center

assert(shift_ratio <= 1);
assert(shift_ratio >= 0);

is_single = isa(image(1),'uint8');

image_size = size(image); 
image_size = image_size(1:2); 

%Find the center of mass of the mask
center_of_mass_mask = findMaskCenterOfMass(image_mask);

%Find the center of the image
image_center = round(image_size/2);

%Find the desired location of the center of mass of the image
desired_location = round((1-shift_ratio)*center_of_mass_mask+shift_ratio*image_center); 

%Find the shift amount
shift_amount = desired_location-center_of_mass_mask; 

%Find the new image bounds
new_image_bounds = [[1;1], image_size.'] + shift_amount.';
new_image_bounds(new_image_bounds < 1) = 1; 

if(new_image_bounds(1,2) > image_size(1))
    new_image_bounds(1,2) = image_size(1);
end

if(new_image_bounds(2,2) > image_size(2))
    new_image_bounds(2,2) = image_size(2);
end

%Find the old image bounds
old_image_bounds = [[1;1], image_size.'] - shift_amount.';
old_image_bounds(old_image_bounds < 1) = 1; 

if(old_image_bounds(1,2) > image_size(1))
    old_image_bounds(1,2) = image_size(1);
end

if(old_image_bounds(2,2) > image_size(2))
    old_image_bounds(2,2) = image_size(2);
end

%Shift the image
image_out = zeros(size(image)); 
image_out(new_image_bounds(1,1):new_image_bounds(1,2),...
    new_image_bounds(2,1):new_image_bounds(2,2),:) = ...
    image(old_image_bounds(1,1):old_image_bounds(1,2),...
    old_image_bounds(2,1):old_image_bounds(2,2),:);

shifted_mask = zeros(size(image_mask));
shifted_mask(new_image_bounds(1,1):new_image_bounds(1,2),...
    new_image_bounds(2,1):new_image_bounds(2,2),:) = ...
    image_mask(old_image_bounds(1,1):old_image_bounds(1,2),...
    old_image_bounds(2,1):old_image_bounds(2,2),:);

if(is_single)
    image_out = uint8(image_out); 
end
end

function [ center_of_mass ] = findMaskCenterOfMass(image_mask)
[row_mask, col_mask] = find(image_mask); 

center_of_mass = [mean(row_mask), mean(col_mask)];
center_of_mass = round(center_of_mass);
end


