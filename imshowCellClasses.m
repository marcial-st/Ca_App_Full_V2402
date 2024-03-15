 function [classImg,NCLS] = imshowCellClasses(NCLS,classImgBase,CELLS,imgLabel)
%% Create Masks per Class %%
[~,jClass]=size(NCLS);
maskClassI = false(size(classImgBase));
i_coffline = 1;
cent_offline = [];
for jMask=1:1:jClass    
    [~,jCV] = size(NCLS(jMask).cVector);
    NCLS(jMask).mask = maskClassI;    
%     NCLS(jMask).color = round(255*rand(1,3));
%     NCLS(jMask).label = strcat("Class",num2str(jMask)); 
    for jMClass =1:1:jCV
        cent = CELLS(NCLS(jMask).cVector(jMClass)).xy_hist(1,:);
%         cent = CELLS(OnCellIdx(NCLS(jMask).cVector(jMClass))).xy_hist(1,:);
        ilabel = imgLabel(cent(2),cent(1));
        if ilabel ~= 0
            NCLS(jMask).mask = logical(NCLS(jMask).mask + (imgLabel == ilabel));
        end
        if CELLS(NCLS(jMask).cVector(jMClass)).status=="OFFLINE"
            cent_offline(i_coffline,1) = cent(1);
            cent_offline(i_coffline,2) = cent(2);
            i_coffline = i_coffline+1;
        end
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
if ~isempty(cent_offline)
    classImg = insertMarker(classImg,cent_offline,'star','color',[0 255 0],'size',5);
end
imshow(classImg)
for jMask=1:1:jClass    
    [~,jCV] = size(NCLS(jMask).cVector);   
    for jMClass =1:1:jCV
        cent = CELLS(NCLS(jMask).cVector(jMClass)).xy_hist(1,:);
        ilabel = imgLabel(cent(2),cent(1));
        % text(cent(1),cent(2),num2str(ilabel),'FontSize',10,'FontWeight','bold','Color',[0,0,0]) % To set global CELLS label 
        text(cent(1),cent(2),num2str(NCLS(jMask).cVector(jMClass)),'FontSize',10,'FontWeight','bold','Color',[0,0,0]) % to set cells_online label
    end
end
colors = reshape([NCLS.color],3,[])';
imlegend(colors/255,[NCLS.label])