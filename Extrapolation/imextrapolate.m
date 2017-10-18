clear variables;
% constants
filename = 'Car_Ad';
ext = '.jpg';

file = sprintf('%s%s',filename, ext);
% imrgb = imresize(imread(file),0.3);
imrgb = imread(file);
sz = max(size(imrgb)) + 100;
imrgb = s_canvasSize(imrgb, [sz, sz],imrgb(1,1,:));
imlab = rgb2lab(imrgb);

options.rotations=[];
options.psz=7;
options.contDistFromBoundary = 3;
options.extrapPixelsSize = 10;

[ ex_rgb,ex_lab] = s_imextrapolate( imrgb, imlab, options );

figure(111);subplot(1,2,1);imshow(ex_rgb);title('after');
figure(111);subplot(1,2,2);imshow(imrgb);title('before');

