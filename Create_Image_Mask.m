%INPUT_IMAGE = 'Car_Ad.jpg'; 
INPUT_IMAGE = 'Perfume.jpg';

image_path = ['../Ad_Images/Images/', INPUT_IMAGE];
in_image = imread(image_path);

figure, imshow(in_image);
h = imfreehand; 
object_mask = createMask(h);
wait(h);
imshow(object_mask);

[~,im_name,ext] = fileparts(image_path);
imwrite(object_mask, ['../Ad_Images/Masks/', im_name, '_Mask', ext]);