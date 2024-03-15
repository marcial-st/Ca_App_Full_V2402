% Image Binarization&Clean
function binImg = bin_and_clean_segnet2(img,segnet_rough)
img=medfilt2(img,[7 7]);
img=adapthisteq(img);
img=adapthisteq(img);


C = semanticseg(img,segnet_rough.net);
binImg = C=='cell';

se = strel('disk',3);
binImg = imopen(binImg,se);
binImg = imfill(binImg,'holes');
se = strel('disk',1);
binImg = imdilate(binImg,se);