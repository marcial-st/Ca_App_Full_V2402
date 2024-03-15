% Image Binarization&Clean
function binImg =bin_and_clean(img,w)
% t=t/10;
% binImg = bradley(img,[w w],t);
% binImg = medfilt2(binImg,[floor(w/4) floor(w/4)]);
padding='replicate';
 img=double(img);
 imgs=img.*img;
meani = averagefilter(img, [w w], padding);
meani2 = averagefilter(imgs, [w w], padding);
va2=meani2-meani.*meani;
stdi2=sqrt(va2);
binImg= true(size(img));
binImg(img < (meani+stdi2)) = 0;
binImg =medfilt2(binImg ,[floor(w/4) floor(w/4)]);
se = strel('disk',3);
binImg = imopen(binImg,se);
binImg = imfill(binImg,'holes');
se = strel('disk',1);
binImg = imdilate(binImg,se);



% Clean image
% se=strel('disk',1,4);
% binImg=imerode(binImg,se);
% binImg=imopen(binImg,se);
% binImg=imerode(binImg,se);
% binImg=imopen(binImg,se);
% % ---
% [L,~]=bwlabel(binImg);
% v=regionprops(L,'Area','Eccentricity');
% aux=struct2cell(v);
% aux=transpose(cell2mat(aux));
% A=aux(:,1);
% B=aux(:,2);
% am=mean(A)
% amin=min(A)
% amax=max(A)
% as=std(A)
% em=mean(B)
% es=std(B)
% PP=find((0.25*am>A)&(A<4*am));
% IB=zeros(size(L));  
% for i=1:1:length(PP)
%    IB=(L==PP(i))+IB;
% end
% binImg = IB;