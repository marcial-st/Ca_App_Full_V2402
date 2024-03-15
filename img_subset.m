function [im340,bin_340,ratio]=img_subset(D340,D380,limInf,limSup,imDim,w,segnet_medium,segnet_rough)
im340=uint8(zeros(imDim(1),imDim(2),limSup-limInf+1));
im380=im340;
bin_340= false(imDim(1),imDim(2),limSup-limInf+1);
ratio=double(zeros(imDim(1),imDim(2),limSup-limInf+1));
h3 = waitbar(0,'Loading image subset & binarization...');
for iInt=0:1:limSup-limInf
    im340(:,:,iInt+1)=imread(strcat(D340(limInf+iInt).folder,'\',D340(limInf+iInt).name));
    im380(:,:,iInt+1)=imread(strcat(D380(limInf+iInt).folder,'\',D380(limInf+iInt).name));
    ratio(:,:,iInt+1)=double(im340(:,:,iInt+1))./double(im380(:,:,iInt+1));
    if (w == 999)
        bin_340(:,:,iInt+1)=bin_and_clean_segnet1(im340(:,:,iInt+1),segnet_medium);
    elseif (w == 998)
        bin_340(:,:,iInt+1)=bin_and_clean_segnet2(im340(:,:,iInt+1),segnet_rough);
    else
        bin_340(:,:,iInt+1)=bin_and_clean(im340(:,:,iInt+1),w);
    end       
    waitbar(iInt/(limSup-limInf))
end
uint8(im340);
logical(bin_340);
close(h3)
