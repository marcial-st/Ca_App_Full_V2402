%==========Clearing================
close all;
clear;
clc;
%==========Paths=========================
basePath{1}='C:\Users\JMartinez\';
idxPath=1;
DB='191017e3';
%-------===============================
%--C:\Users\JMartinez\Google Drive\CaExp_2020\191017 OZDF (completo)\191017e1
switch DB
        
    case '191017e1'
        % %============DS_1 Zheng===============
        s_path=[ basePath{idxPath} 'Google Drive\' 'CaExp_2020\' DB '\340\'];
        wPath=[s_path 'Seg'];
        ftype='png';
        nlevels=3;
        %         mProm=ceil(120*120/1000);
        
    case '191017e3'
        % %============DS_2 CellaVision Zheng===============
        s_path=[ basePath{idxPath} 'Google Drive\' 'CaExp_2020\' DB '\340\'];
        wPath=[s_path 'Seg'];
        ftype='png';
        nlevels=3;
        %         mProm=300*300/1000;
    case '190604e3'
        s_path=[ basePath{idxPath} 'Google Drive\' 'CaExp_2020\190604 OZDF (completo)\' DB '\380\'];
        wPath=[s_path 'Seg'];
        ftype='png';
        nlevels=2;
    case '190607e5'
        s_path=[ basePath{idxPath} 'Google Drive\' 'CaExp_2020\190607 OZDF (completo)\' DB '\380\'];
        wPath=[s_path 'Seg'];
        ftype='png';
        nlevels=3;
end
%========Begin code=======
% mPath=w_path{1};
dfiles=dir([s_path '*.' ftype] ); %1 to 999
% dfiles=dir([s_path '*-1???.' ftype] ); %1000 to 9999
lastIm=length(dfiles);

%==========Variable defaults================
chLevel=0:255;
intervals=7.0;
winSz=9;

for mImg=246:250%lastIm
    ofname=dfiles(mImg).name;
    oCells1=imread([s_path  ofname]);
    oCells=ordfilt2(oCells1,floor(winSz*winSz/1.34),ones(winSz,winSz));
%     mxVal=max(max(oCells));
%     inFactor=floor(255/mxVal);
    %     oCells=inFactor*oCells;
    %     scols=sum(oCells);
    %     figure;
    %     plot(scols);
    %     oCells=adapthisteq(oCells1,'NumTiles',[15 21]);
    %     figure;
    %     imshow(oCells*3);
    %     imshowpair(oCells,oCells1,'montage');
    % figure;
    % imhist(oCells);
    %
    %     imAdj=imadjust(oCells);
    %     figure;
    %     imshow(imAdj);
    nrows=size(oCells,1);
    ncols=size(oCells,2);
    bands=[0,floor(ncols*0.10), floor(ncols*0.825), ncols];
    b0=bands(2)+1;
    b1=bands(3);
    %     imRes=uint8(zeros(nrows,ncols,3));
    imRes=uint8(zeros(nrows,b1-b0+1,3));
    
    %     imRes= imRes+250*uint8(oCells > 39) +200*uint8(oCells <= 39 & oCells > 36)+ 150* uint8(oCells <= 36 & oCells > 30);
    %         figure('Name', 'CaFlow');
    %         imshow(oCells > 37);
    %     oCells=imAdj;
    % figure;
    his=imhist(oCells(:,b0:b1));
    %     plot(his);
    %                 plot(0:50,his(1:51));
    [vMx,idmx]=max(his);
    %             TbPks=islocalmax(his,'MinSeparation',3,'MinProminence',100);
    %             Tbmin=islocalmin(his,'MinSeparation',3,'MinProminence',100);
    %             bPks=chLevel(TbPks);
    %             bMins=chLevel(Tbmin);
    %           centroids=bPks(end-6:end)';
    hlast=find (his(idmx:end) ==0,1,'first')+idmx-1;
    delta= floor(double(hlast- idmx)/intervals);
    if delta <1
        delta=1;
    end
    dIntensity=floor(250.0/intervals);
    centroids(1,1)=10;
    for ii=2:intervals+4
        centroids(ii,1)=idmx + (ii-6)*delta;
    end
    %     ab = reshape(double(oCells),nrows*ncols,1);
    ab = reshape(double(oCells(:,b0:b1)),nrows*(b1-b0+1),1);
    
    nClus=length(centroids);
    [cluster_idx, C]=kmeans(ab,nClus, 'Distance','cityblock', ...
        'Start',centroids);
    [Cs, idC]=sort(C);
    %     pixel_labels = reshape(cluster_idx,nrows,ncols);
    pixel_labels = reshape(cluster_idx,nrows,(b1-b0+1));
    rFac=floor(255/nClus);
    for k=nClus-(nlevels-1):nClus
        imRes(:,:,1)=imRes(:,:,1) + uint8((pixel_labels==idC(k)))*k*rFac;
        %             if k==2
        %                 imW=imbinarize(imRes(:,:,channel));
        %                 imF=imfill(imW,'holes');
        %                 figure('Name',color{channel});
        %                 imshow(imF);
        %             end
    end
    gFac=floor(255/(nClus-nlevels));
    
    for k=5:nClus-(nlevels)
        imRes(:,:,2)=imRes(:,:,2) + uint8((pixel_labels==idC(k)))*k*gFac;
    end
    bFac=floor(255/4);
    
    for k=2:4
        imRes(:,:,3)=imRes(:,:,3) + uint8((pixel_labels==idC(k)))*k*bFac;
    end
    figure('Name', [ofname '_' DB '_' 'CaFlow']);
    %     imshow(imRes);
    imshowpair(oCells(:,b0:b1),imRes,'montage');
end