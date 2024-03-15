function [im340,bin_340,ratio]=img_subset_plist(D340,D380,limInf,limSup,imDim,process_list,segnet_medium,segnet_rough)
im340=uint8(zeros(imDim(1),imDim(2),limSup-limInf+1));
im380=im340;
bin_340= false(imDim(1),imDim(2),limSup-limInf+1);
ratio=double(zeros(imDim(1),imDim(2),limSup-limInf+1));
h3 = waitbar(0,'Loading image subset & binarization...');
for iInt=0:1:limSup-limInf
    im340(:,:,iInt+1)=imread(strcat(D340(limInf+iInt).folder,'\',D340(limInf+iInt).name));
    im380(:,:,iInt+1)=imread(strcat(D380(limInf+iInt).folder,'\',D380(limInf+iInt).name));
    ratio(:,:,iInt+1)=double(im340(:,:,iInt+1))./double(im380(:,:,iInt+1));

    current_img_p = im340(:,:,iInt+1);
    
    idx = max(size(process_list));    
    for i =1:1:idx
        selection = string(process_list(i).name);
        switch selection
            case 'MedianFilter'
                current_img_p = medfilt2(current_img_p,process_list(i).params);
                %%disp(selection)
            case 'Contrast Eq'
                current_img_p = histeq(current_img_p);
                %disp(selection)
            case 'Adaptive Contrast Eq'
                current_img_p = adapthisteq(current_img_p);
                %disp(selection)
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
                %disp(selection)
            case 'SegNet Medium'
                C = semanticseg(current_img_p,segnet_medium.net);
                current_img_p = C=='cell';
                current_img_p(1:10,:) = 0;
                current_img_p(:,1:10) = 0;
                current_img_p(end-10:end,:) = 0;
                current_img_p(:,end-10:end) = 0;
                logical(current_img_p);                             
                %disp(selection)
            case 'SegNet Rough'
                C = semanticseg(current_img_p,segnet_rough.net);
                current_img_p = C=='cell';
                current_img_p(1:10,:) = 0;
                current_img_p(:,1:10) = 0;
                current_img_p(end-10:end,:) = 0;
                current_img_p(:,end-10:end) = 0;
                logical(current_img_p);
                %disp(selection)                                
            case 'Erode'
                se = strel('disk',process_list(i).params);
                current_img_p = imerode(current_img_p,se);
                %disp(selection)
            case 'Dilate'
                se = strel('disk',process_list(i).params);
                current_img_p = imdilate(current_img_p,se);
                %disp(selection)
            case 'OpenFilt'
                se = strel('disk',process_list(i).params);
                current_img_p = imopen(current_img_p,se);
                %disp(selection)
            case 'CloseFilt'
                se = strel('disk',process_list(i).params);
                current_img_p = imclose(current_img_p,se);
                %disp(selection)
            case 'FillHoles'
                current_img_p = imfill(current_img_p,'holes');
                %disp(selection)
            case 'AreaThreshold'
                if islogical(current_img_p)
                    if isempty(process_list(i).params)
                        stats = regionprops(current_img_p,'Area');
                        figure;histogram([stats.Area]);title('Area histogram');xlabel('Area');ylabel('Frequency');
                        range = processDefAreaFilter;
                        process_list(i).params = range;                    
                    end                    
                    current_img_p = bwareafilt(current_img_p,process_list(i).params);                    
                else
                    warndlg('Area filter only for binary images','Warning');
                end
                %disp(selection)
            otherwise
                %disp('holis')
        end
    end
    current_img_p(1:10,:) = 0;
    current_img_p(:,1:10) = 0;
    current_img_p(end-10:end,:) = 0;
    current_img_p(:,end-10:end) = 0;
    
    bin_340(:,:,iInt+1) = current_img_p;
    
    waitbar(iInt/(limSup-limInf))
end
uint8(im340);
logical(bin_340);
close(h3)