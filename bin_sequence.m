%%% Thesis
% Adaptive Binarization Function
% Bradley Technique
% Median filtering with 7x7 window is needed before this function
function [b_seq,m,n,k,imb]=bin_sequence(in_seq,w)
[m,n,k]=size(in_seq);
b_seq=zeros(m,n,k);
l=floor(w/2);
h3 = waitbar(0,'Adaptive binarization progress...');
for i=1:1:k
    im340f=adapthisteq(medfilt2(in_seq(:,:,i))); 
    im340f=im2double(medfilt2(im340f,[7,7]));
    b_seq(:,:,i)=bin_and_clean(im340f,w);
%     b_seq(:,:,i)=ac_bradley(imb);
    waitbar(i / k)
end
close(h3)
% function [b_seq,m,n,k,imb]=bin_sequence(in_seq,t,w)
% [m,n,k]=size(in_seq);
% b_seq=zeros(m,n,k);
% l=floor(w/2);
% h3 = waitbar(0,'Adaptive binarization progress...');
% for i=1:1:k
%     im340f=adapthisteq(medfilt2(in_seq(:,:,i))); 
%     im340f=im2double(medfilt2(im340f,[7,7]));
% %     imb=not(bradley(im340f,t,w));
%     imb=not(bradley(im340f,[w w],t));    
%     [counts,x]=imhist(imb(l:m-l,l:n-l));
%     i_pad=find(counts==max(counts));
%     pad=logical(x(i_pad));
%     imb(1:l+1,:)=pad;
%     imb(m-l:m,:)=pad;
%     imb(:,1:l+1)=pad;
%     imb(:,n-l:n)=pad;
%     b_seq(:,:,i)=ac_bradley(imb);
%     waitbar(i / k)
% end
% close(h3)