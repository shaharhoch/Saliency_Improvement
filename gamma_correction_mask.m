function [ out_img ] = gamma_correction_mask( in_image, gamma_mask )
in_img_double = im2double(in_image);

out_img = zeros(size(in_image));

for ind=1:size(in_image, 3)
    out_img(:,:,ind) = in_img_double(:,:,ind) .^ gamma_mask;
end

if(isa(in_image,'single'))
    out_img = im2single(out_img);
end
end

