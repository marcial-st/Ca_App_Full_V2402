% Image Binarization&Clean
function binImg = bin_and_clean_segnet1(img,segnet_medium)
C = semanticseg(img,segnet_medium.net);
binImg = C=='cell';
binImg(1:10,:) = 0;
binImg(:,1:10) = 0;
binImg(end-10:end,:) = 0;
binImg(:,end-10:end) = 0;
se = strel('disk',4);
binImg = imopen(binImg,se);
binImg = imfill(binImg,'holes');
% se = strel('disk',1);
% binImg = imdilate(binImg,se);