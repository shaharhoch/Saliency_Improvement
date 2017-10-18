function [ mask ] = ImageBoardersToMask( image, image_boarders )
% Gets image boarders and returns an image mask
mask = zeros(size(image,1), size(image,2));
mask(image_boarders(1,1):image_boarders(1,2),image_boarders(2,1):image_boarders(2,2)) = 1;

end

