clear

%% parameters
testImageName = 'Car_Ad'; % cow or bungee or man
psz = 9; % patch size ( to be inpainted in inpainting)
testImagePath = '../../Ad_Images/Images/';%'~/Documents/MATLAB/AutoShared/testimages/Petter_Strandmark/';
testImageSource = fullfile(testImagePath,testImageName);
origImg = imread([testImageSource,'.jpg']);

mask_path = '../../Ad_Images/Masks/';
mask = imread([mask_path,testImageName,'_Mask.jpg']);
mask(mask==255) = 1;
mask = logical(mask);

tic
%output is fix!
[inpaintedImg,c,d,fillingMovie] = inpainting(origImg,mask,psz);
toc
maskedImg = repmat(uint8(~mask),[1,1,3]).*origImg;

figure(1),imshow(uint8(origImg)),title('Original Image')
figure(2),imshow(uint8(maskedImg)),title('Masked Image')
figure(3),imshow(uint8(inpaintedImg)),title('Inpainted Image')
folderName = ['myresults/',datestr(now,'yymmdd-HHMMSS'),'_',testImageName];
mkdir(folderName)

imwrite(uint8(origImg),fullfile(folderName,'origImg.bmp'),'BMP');
imwrite(uint8(maskedImg),fullfile(folderName,'maskedImg.bmp'),'BMP');
imwrite(uint8(inpaintedImg),fullfile(folderName,'inpaintedImg.bmp'),'BMP');
