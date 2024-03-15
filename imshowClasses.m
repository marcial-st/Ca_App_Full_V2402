function imshowClasses(NCLS,classImgBase,CELLS,OnCellIdx,imgLabel)
%% Create Masks per Class %%
[~,jClass]=size(NCLS);
maskClassI = false(size(classImgBase));
for jMask=1:1:jClass
    [~,jCV] = size(NCLS(jMask).cVector);
    NCLS(jMask).mask = maskClassI;    
    NCLS(jMask).color = round(255*rand(1,3));
    NCLS(jMask).label = strcat("Class ",num2str(jMask)); 
    for jMClass =1:1:jCV
        cent = CELLS(OnCellIdx(NCLS(jMask).cVector(jMClass))).xy_hist(1,:);
        ilabel = imgLabel(cent(2),cent(1));
        NCLS(jMask).mask = logical(NCLS(jMask).mask + (imgLabel == ilabel));
    end
    classImgBase = classImgBase.*(uint8(~NCLS(jMask).mask));
end
%%% Show spatial location per class
classImg = uint8(zeros([size(classImgBase) 3]));
classImg(:,:,1)=classImgBase;
classImg(:,:,2)=classImgBase;
classImg(:,:,3)=classImgBase;
for jMask=1:1:jClass
    classImg(:,:,1)=classImg(:,:,1)+(NCLS(jMask).color(1,1).*uint8(NCLS(jMask).mask));
    classImg(:,:,2)=classImg(:,:,2)+(NCLS(jMask).color(1,2).*uint8(NCLS(jMask).mask));
    classImg(:,:,3)=classImg(:,:,3)+(NCLS(jMask).color(1,3).*uint8(NCLS(jMask).mask));
end
figure;imshow(classImg)
colors = reshape([NCLS.color],3,[])';
imlegend(colors/255,[NCLS.label])