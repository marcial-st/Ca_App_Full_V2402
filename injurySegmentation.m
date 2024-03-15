function [bin_img,inj_img] = injurySegmentation(img)
w=40;
[i_img,j_img] = size(img);
imed=medfilt2(img,[w w]);
imeda2=adapthisteq(imed);
% figure;
% % hold(f1,'on')
% hold on
space=[1:1:j_img];
imgaux=imeda2;
injMask=false(size(img));
f = waitbar(0,'Thanks for waiting...');
% pause(0.01)
bar_message = "Injury segmentation in process...";
for k=1:1:i_img-2
    imgaux=imeda2;
    imgaux(k,:)=1;
    vector=smooth(double(imeda2(k,:)));
    [TF,P] = islocalmin(vector,'MinProminence',3);
    % subplot(1,2,1)    
    % yyaxis left
    % plot(space,vector,space(TF),vector(TF),'r*')
    % yyaxis right
    % plot(space,P,'g','DisplayName','Prominence')
    % legend show
    spAux=space(TF)';
    pAux=floor(P(TF)/2);    
    if ~isempty(spAux)
        [iMask,~]=size(spAux);
        for ifm=1:1:iMask
            injMask(k,(spAux(ifm)-pAux(ifm)):(spAux(ifm)+pAux(ifm)))=true;
        end    
    end
    % subplot(1,2,2)
    % number=num2str(k);   
    % imshow(imgaux)
    % title(strcat('Image row: ',number))
    % pause(0.05)
    % cla
    waitbar(k/(i_img-2),f,bar_message);
end

bin_img=injMask;
se=strel('disk',2,6);
bin_img=imerode(bin_img,se);
bin_img=imerode(bin_img,se);
bin_img=imdilate(bin_img,se);
bin_img=imdilate(bin_img,se);

inj_img=uint8(zeros(i_img,j_img,3));
inj_img(:,:,1)=img+(255*uint8(bin_img));
inj_img(:,:,2)=img;
inj_img(:,:,3)=img;
close(f)

% results=figure;
% subplot(1,2,1)
% imshow(injMask)
% subplot(1,2,2)
% imshow(inj_img)