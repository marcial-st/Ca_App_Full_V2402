function current_img_p = applyProcessList(process_list,current_img_p)

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
            disp('holis')
    end
end