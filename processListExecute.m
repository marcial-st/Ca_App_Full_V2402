function [current_img_p] = processListExecute(current_img,process_list)
%PROCESSLISTEXECUTE Summary of this function goes here
%   Detailed explanation goes here
% global process_list
% global current_img
% global current_img_p
% global segnet_medium
% global segnet_rough
segnet_medium = load('SegNet_medium_v200714.mat');
segnet_rough = load('SegNet_rough_v200527.mat');
% process_dir = {'--- Enhancement ---',...
%                'MedianFilter','Contrast Eq','Adaptive Contrast Eq',...
%                '--- Binarization ---',...
%                'Adaptative Bin','SegNet Medium','SegNet Rough',...
%                '--- Post Bin ---',...
%                'Erode','Dilate','OpenFilt','CloseFilt','AreaThreshold'};
current_img_p = current_img;
idx = max(size(process_list));
if idx ~= 0    
    for i =1:1:idx
        selection = string(process_list(i).name);
        switch selection
            case 'MedianFilter'
                current_img_p = medfilt2(current_img_p,process_list(i).params);
                disp(selection)
            case 'Contrast Eq'
                current_img_p = histeq(current_img_p);
                disp(selection)
            case 'Adaptive Contrast Eq'
                current_img_p = adapthisteq(current_img_p);
                disp(selection)
            case 'Adaptative Bin'
                w = process_list(i).params;
                padding='replicate';
                img=double(current_img_p);
                imgs=img.*img;
                meani = averagefilter(img, [w w], padding);
                meani2 = averagefilter(imgs, [w w], padding);
                va2=meani2-meani.*meani;
                stdi2=sqrt(va2);
                current_img_p= true(size(img));
                current_img_p(img < (meani+stdi2)) = 0;
                logical(current_img_p);
                disp(selection)
            case 'SegNet Medium'
                C = semanticseg(current_img_p,segnet_medium.net);
                current_img_p = C=='cell';
                current_img_p(1:10,:) = 0;
                current_img_p(:,1:10) = 0;
                current_img_p(end-10:end,:) = 0;
                current_img_p(:,end-10:end) = 0;
                logical(current_img_p);                             
                disp(selection)
            case 'SegNet Rough'
                C = semanticseg(current_img_p,segnet_rough.net);
                current_img_p = C=='cell';
                current_img_p(1:10,:) = 0;
                current_img_p(:,1:10) = 0;
                current_img_p(end-10:end,:) = 0;
                current_img_p(:,end-10:end) = 0;
                logical(current_img_p);
                disp(selection)
            case 'Kmeans'
                
                chLevel=0:255;
%                 intervals=7.0;
                intervals = process_list(i).params;
                winSz=9;
                nlevels=3;
    
                nrows=size(current_img_p,1);
                ncols=size(current_img_p,2);
%                 bands=[0,floor(ncols*0.10), floor(ncols*0.825), ncols];
%                 b0=bands(2)+1;
%                 b1=bands(3);
%                 imRes=uint8(zeros(nrows,b1-b0+1,3));
                imRes=uint8(zeros(nrows,ncols,3));
    
                his=imhist(current_img_p(:,:));
%                 his=imhist(current_img_p(:,b0:b1));
                [vMx,idmx]=max(his);
                hlast=find(his(idmx:end) ==0,1,'first')+idmx-1;
                delta= floor(double(hlast-idmx)/intervals);
                if delta <1
                    delta = 1;
                end
                dIntensity=floor(250.0/intervals);
                centroids(1,1)=10;
                for ii=2:intervals+4
                    centroids(ii,1)=idmx + (ii-6)*delta;
                end
                ab = reshape(double(current_img_p),nrows*ncols,1);
%                 ab = reshape(double(current_img_p(:,b0:b1)),nrows*(b1-b0+1),1);

                nClus=length(centroids);
                [cluster_idx, C]=kmeans(ab,nClus, 'Distance','cityblock', ...
                    'Start',centroids);
                [Cs, idC]=sort(C);
%                 pixel_labels = reshape(cluster_idx,nrows,(b1-b0+1));
                pixel_labels = reshape(cluster_idx,nrows,ncols);
                rFac=floor(255/nClus);
                for k=nClus-(nlevels-1):nClus
                    imRes(:,:,1)=imRes(:,:,1) + uint8((pixel_labels==idC(k)))*k*rFac;
                end
                gFac=floor(255/(nClus-nlevels));
                for k=5:nClus-(nlevels)
                    imRes(:,:,2)=imRes(:,:,2) + uint8((pixel_labels==idC(k)))*k*gFac;
                end
                bFac=floor(255/4);
                for k=2:4
                    imRes(:,:,3)=imRes(:,:,3) + uint8((pixel_labels==idC(k)))*k*bFac;
                end
                current_img_p = imRes(:,:,1);
                disp(selection)
            case 'Erode'
                se = strel('disk',process_list(i).params);
                current_img_p = imerode(current_img_p,se);
                disp(selection)
            case 'Dilate'
                se = strel('disk',process_list(i).params);
                current_img_p = imdilate(current_img_p,se);
                disp(selection)
            case 'OpenFilt'
                se = strel('disk',process_list(i).params);
                current_img_p = imopen(current_img_p,se);
                disp(selection)
            case 'CloseFilt'
                se = strel('disk',process_list(i).params);
                current_img_p = imclose(current_img_p,se);
                disp(selection)
            case 'FillHoles'
                current_img_p = imfill(current_img_p,'holes');
                disp(selection)
            case 'AreaThreshold'
                if islogical(current_img_p)
                    if isempty(process_list(i).params)
                        stats = regionprops(current_img_p,'Area');
                        figure;histogram([stats.Area]);title('Area histogram');xlabel('Area');ylabel('Frequency');
                        range = processDefAreaFilter;
                        process_list(i).params = range;                    
                    end                    
                    current_img_p = bwareafilt(current_img_p,process_list(i).params);
%                     processListUpdate(handles)
                else
                    warndlg('Area filter only for binary images','Warning');
                end
                disp(selection)
            otherwise
                disp('holis')
        end
    end
    if islogical(current_img_p)
        current_img_p(1:10,:) = 0;
        current_img_p(:,1:10) = 0;
        current_img_p(end-10:end,:) = 0;
        current_img_p(:,end-10:end) = 0;
        [m,n]=size(current_img_p);
        imgrgb=uint8(zeros(m,n,3));
        imgrgb(:,:,1)=uint8(current_img); imgrgb(:,:,2)=imgrgb(:,:,1); imgrgb(:,:,3)=imgrgb(:,:,1)+uint8(current_img_p).*255;
%         figure
%         axes(handles.axes1)%pongo la imagen en la ventana uno
%         imshow(imgrgb) 
    else
%         imshow(current_img_p)
    end    
else
%     imshow(current_img)
end
end

