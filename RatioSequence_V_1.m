%%% Thesis
%%% Ratio Sequence V_1
%%% Load and save 340 nm image sequence
function [im340,im380,ratio,n_ratio]=RatioSequence_V_1

[filename, pathname] = uigetfile('*.tif','Select the 340nm first image');
if filename~=0
    D = dir(strcat([pathname, '*.tif']));
    nf = length(D(not([D.isdir])));
    guion=strfind(filename,'-');
    punto=strfind(filename,'.');
    ni=str2num(filename(punto-3:punto-1));
    h1 = waitbar(0,'Loading...');
    for i=ni:1:ni+nf-1
        nim=strcat(pathname,filename(1:guion),sprintf('%03d',i),filename(punto:length(filename)));
        im340(:,:,i-ni+1)=imread(nim);
        waitbar(i / nf)
    end
    close(h1);
%%% Load and save 380 nm image sequence
    [filename, pathname] = uigetfile('*.tif','Select the 380nm first image');
    if filename ==0
        warndlg('Missing one image sequence','Warning !')
        im340=0;
        im380=0;
        ratio=0;
        n_ratio=0;
        return
    else
        D = dir(strcat([pathname, '*.tif']));
        nf = length(D(not([D.isdir])));
        guion=strfind(filename,'-');
        punto=strfind(filename,'.');
        ni=str2num(filename(punto-3:punto-1));
        h2 = waitbar(0,'Loading...');
        for i=ni:1:ni+nf-1
            nim=strcat(pathname,filename(1:guion),sprintf('%03d',i),filename(punto:length(filename)));
            im380(:,:,i-ni+1)=imread(nim);
            waitbar(i / nf)
        end
        close(h2);
    end   
%%% Check array dimensions match
    s340=size(im340);
    s380=size(im380);
    if s340~=s380
     warndlg('340nm and 380nm size missmatch.','Warning !')
      flag_c=0;
    else
        flag_c=1;
    end;
%% Makes Ratio image sequence
    ratio=zeros(s340);
    if 1==flag_c
      for i=1:1:s340(length(s340))    
          ratio(:,:,i)=(double(im340(:,:,i))./double(im380(:,:,i)));
      end
      max_r=max(max(max(ratio(:,:,:))));
      n_ratio=uint8(255*((ratio./max_r)));
      max_n=max(max(max(n_ratio(:,:,:))));
      min_n=min(min(min(n_ratio(:,:,:))));
     n_ratio=(255/(255-min_n))*(n_ratio-min_n);
    end
%% Contrast Enhacenment in Ratio sequence
else
    im340=0;
    im380=0;
    ratio=0;
    n_ratio=0;
end